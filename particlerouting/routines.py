# -*- coding: utf-8 -*-
"""
Scrips to simplify process of routing and visualizing tracer particles.

Project Homepage: https://github.com/
"""

from .particle_track import Particle
from math import floor, sqrt, pi
import numpy as np
from random import shuffle
import matplotlib
from matplotlib import pyplot as plt
from scipy import ndimage
import sys, os, re, string
from netCDF4 import Dataset
import time as time_lib
from scipy.sparse import lil_matrix, csc_matrix, hstack
import logging
import time



def steady_plots(params,num_iter,folder_name):
    '''
    Function to automate plotting of particle movement over a set of steady
    fields (time-invariant)

    Inputs :
                params : class of parameter values for the particles
                num_iter : number of iterations to move particles over
                folder_name : string of folder name to put outputs in

    Outputs :
                Script saves result of each iteration to a folder
                with both the figure for each iteration as a png and the data
                with the particle start and end locations
    '''

    # define the particle
    particle = Particle(params)

    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except:
        print('Directories already exist')

    # iterate and save results
    for i in range(0,num_iter):
        # do particle iteration
        if i == 0:
            start_inds, end_inds, travel_times = particle.run_iteration()
            beg_inds = start_inds # for 1st timestep
            xinds=[];yinds=[];
            for j in range(0,len(end_inds)):
                xinds.append(end_inds[j][0])
                yinds.append(end_inds[j][1])
        else:
            beg_inds, end_inds, travel_times = particle.run_iteration(start_xindices=xinds,start_yindices=yinds,start_times=travel_times)
            xinds = []; yinds = [];
            for j in range(0,len(end_inds)):
                xinds.append(end_inds[j][0])
                yinds.append(end_inds[j][1])

        # do plot and saving
        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,len(start_inds)):
            plt.scatter(start_inds[k][1],start_inds[k][0],c='b',s=0.75)
            plt.scatter(end_inds[k][1],end_inds[k][0],c='r',s=0.75)
        plt.imshow(params.depth)
        plt.title('Depth')
        plt.colorbar()
        plt.savefig(os.getcwd() + '/' + folder_name + '/figs/output'+str(i)+'.png')
        plt.close()

        # save data
        np.savez(
                os.getcwd() + '/' + folder_name + '/data/data'+str(i)+'.npz',
                beg_inds=beg_inds,
                end_inds=end_inds,
                travel_times=travel_times
                )


def unsteady_plots(params, num_steps, timestep,
                       output_base, output_type,
                       first_output, last_output,
                       folder_name):
    '''
    Function to automate plotting of particle movement in an unsteady flow
    field (time-varying). Each particle travels for the length of the timestep
    before the next output is loaded and the flow field is updated.

    Inputs :
                params : class of particle parameter values
                num_steps : number of model timesteps being covered
                timestep : model timestep duration (seconds)
                output_base : filepath where hydrodynamic outputs are
                output_type : convention output files have been saved as
                              (limited built-in support, may require modification)
                first_output : number to append to output_base for first file
                last_output : number to append to output_base for last file
                folder_name : string of the desired output folder name

    Outputs :
                Script saves result of each iteration to a folder
                with both the figure for each iteration as a png and the data
                with the particle start and end locations
    '''

    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except:
        print('Directories already exist')

    for i in range(first_output,last_output+1):
        # load depth, qx, qy for this timestep
        # making assumption that other variables are constant between output files !!!!
        if output_type == 'csv':
            depth = np.loadtxt(output_base+'/depth'+str(i)+'.csv', delimiter=',')
            qx = np.loadtxt(output_base+'/qx'+str(i)+'.csv', delimiter=',')
            qy = np.loadtxt(output_base+'/qy'+str(i)+'.csv', delimiter=',')
        elif output_type == 'npy':
            depth = np.load(output_base+'/depth'+str(i)+'.npy')
            qx = np.load(output_base+'/qx'+str(i)+'.npy')
            qy = np.load(output_base+'/qy'+str(i)+'.npy')
        elif output_type == 'npz':
            data = np.load(output_base+'/data'+str(i)+'.npz')
            depth = data['depth']
            qx = data['qx']
            qy = data['qy']
        else:
            raise ValueError('Output datatype/structure unsupported, modify the output reading portion of the code')

        # then define the particle class and continue
        particle = Particle(params)

        # reset time
        t = 0.0
        while t < timestep*num_steps:
            if i == first_output and t == 0.0:
                t_old = 0.0
                start_inds, end_inds, travel_times = particle.run_iteration(time_step=timestep)
                xinds=[];yinds=[];
                for j in range(0,len(end_inds)):
                    xinds.append(end_inds[j][0])
                    yinds.append(end_inds[j][1])
            else:
                t_old = np.mean(travel_times)
                beg_ind, end_inds, travel_times = particle.run_iteration(start_xindices=xinds,start_yindices=yinds,start_times=travel_times,time_step=timestep)
                xinds = []; yinds = [];
                for j in range(0,len(end_inds)):
                    xinds.append(end_inds[j][0])
                    yinds.append(end_inds[j][1])

            print('mean travel time: ' + str(np.mean(travel_times)-t_old))
            t = np.mean(travel_times) - t_old

        # make and save plots and data
        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,len(start_inds)):
            plt.scatter(start_inds[k][1],start_inds[k][0],c='b',s=0.75)
            plt.scatter(end_inds[k][1],end_inds[k][0],c='r',s=0.75)
        cbar = plt.colorbar()
        cbar.set_label('Particle Travel Times [s]')
        plt.imshow(params.depth)
        plt.title('Depth')
        cbar2 = plt.colorbar(orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.savefig(os.getcwd() + '/' + folder_name + '/figs/output'+str(i)+'.png')
        plt.close()

        # save data
        np.savez(
                os.getcwd() + '/' + folder_name + '/data/data'+str(i)+'.npz',
                beg_inds=beg_inds,
                end_inds=end_inds,
                travel_times=travel_times
                )


def unsteady_avg_plots(params, avg_timestep,
                       output_base, output_type,
                       first_output, last_output,
                       folder_name):
    '''
    Function to automate plotting of particle movement in an unsteady flow
    field (time-varying). Steps through time as an average of the particle
    travel times.

    Inputs :
                params : class of particle parameter values
                avg_timestep : average timestep duration (estimated using the
                                average particle travel times)
                output_base : filepath where hydrodynamic outputs are
                output_type : convention output files have been saved as
                              (limited built-in support, may require modification)
                first_output : number to append to output_base for first file
                last_output : number to append to output_base for last file
                folder_name : string of the desired output folder name

    Outputs :
                Script saves result of each iteration to a folder
                with both the figure for each iteration as a png and the data
                with the particle start and end locations
    '''

    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except:
        print('Directories already exist')

    for i in range(first_output,last_output+1):
        # load depth, qx, qy for this timestep
        # making assumption that other variables are constant between output files !!!!
        if output_type == 'csv':
            depth = np.loadtxt(output_base+'/depth'+str(i)+'.csv', delimiter=',')
            qx = np.loadtxt(output_base+'/qx'+str(i)+'.csv', delimiter=',')
            qy = np.loadtxt(output_base+'/qy'+str(i)+'.csv', delimiter=',')
        elif output_type == 'npy':
            depth = np.load(output_base+'/depth'+str(i)+'.npy')
            qx = np.load(output_base+'/qx'+str(i)+'.npy')
            qy = np.load(output_base+'/qy'+str(i)+'.npy')
        elif output_type == 'npz':
            data = np.load(output_base+'/data'+str(i)+'.npz')
            depth = data['depth']
            qx = data['qx']
            qy = data['qy']
        else:
            raise ValueError('Output datatype/structure unsupported, modify the output reading portion of the code')

        # then define the particle class and continue
        particle = Particle(params)

        # reset time
        t = 0.0
        while t < avg_timestep:
            if i == first_output and t == 0.0:
                t_old = 0.0
                start_inds, end_inds, travel_times = particle.run_iteration()
                xinds=[];yinds=[];
                for j in range(0,len(end_inds)):
                    xinds.append(end_inds[j][0])
                    yinds.append(end_inds[j][1])
            else:
                t_old = np.mean(travel_times)
                beg_ind, end_inds, travel_times = particle.run_iteration(start_xindices=xinds,start_yindices=yinds,start_times=travel_times)
                xinds = []; yinds = [];
                for j in range(0,len(end_inds)):
                    xinds.append(end_inds[j][0])
                    yinds.append(end_inds[j][1])

            print('mean travel time: ' + str(np.mean(travel_times)-t_old))
            t = np.mean(travel_times) - t_old

        # make and save plots and data
        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,len(start_inds)):
            plt.scatter(start_inds[k][1],start_inds[k][0],c='b',s=0.75)
            plt.scatter(end_inds[k][1],end_inds[k][0],c='r',s=0.75)
        cbar = plt.colorbar()
        cbar.set_label('Particle Travel Times [s]')
        plt.imshow(params.depth)
        plt.title('Depth')
        cbar2 = plt.colorbar(orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.savefig(os.getcwd() + '/' + folder_name + '/figs/output'+str(i)+'.png')
        plt.close()

        # save data
        np.savez(
                os.getcwd() + '/' + folder_name + '/data/data'+str(i)+'.npz',
                beg_inds=beg_inds,
                end_inds=end_inds,
                travel_times=travel_times
                )



def time_plots(params,num_iter,folder_name):
    '''
    Make plots with each particle's travel time visualized.
    Routine assumes a steady flow field, but could be expanded to an unsteady case.

    Inputs :
                params : parameters for the particle
                num_iter : number of iterations to move particles
                folder_name : string of desired output folder name

    Outputs :
                Saves plots and data for each iteration
    '''

    # define the particle
    particle = Particle(params)

    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except:
        print('Directories already exist')

    # iterate and save results
    for i in range(0,num_iter):
        # do particle iteration
        if i == 0:
            start_inds, end_inds, travel_times = particle.run_iteration()
            beg_inds = start_inds # for 1st timestep
            xinds=[];yinds=[];
            for j in range(0,len(end_inds)):
                xinds.append(end_inds[j][0])
                yinds.append(end_inds[j][1])
        else:
            beg_inds, end_inds, travel_times = particle.run_iteration(start_xindices=xinds,start_yindices=yinds,start_times=travel_times)
            xinds = []; yinds = [];
            for j in range(0,len(end_inds)):
                xinds.append(end_inds[j][0])
                yinds.append(end_inds[j][1])

        # do plot and saving
        plt.figure(figsize=(4,4),dpi=200)
        cm = matplotlib.cm.colors.Normalize(vmax=np.max(travel_times), vmin=np.min(travel_times))
        for k in range(0,len(start_inds)):
            plt.scatter(start_inds[k][1],start_inds[k][0],c='b',s=0.75)
            plt.scatter(end_inds[k][1],end_inds[k][0],c=travel_times[k],s=0.75,cmap='coolwarm',norm=cm)
        cbar = plt.colorbar()
        cbar.set_label('Particle Travel Times [s]')
        plt.imshow(params.depth)
        plt.title('Depth')
        cbar2 = plt.colorbar(orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.savefig(os.getcwd() + '/' + folder_name + '/figs/output'+str(i)+'.png')
        plt.close()

        # save data
        np.savez(
                os.getcwd() + '/' + folder_name + '/data/data'+str(i)+'.npz',
                beg_inds=beg_inds,
                end_inds=end_inds,
                travel_times=travel_times
                )



### Function to automate animation of the png outputs
# requires installation of the animation writer 'ffmpeg' which is not part of the
# default PartingRouting installation set of packages
def animate_plots(start_val,end_val,folder_name):
    '''
    Routine to make mp4 animation of the particle routing from png outputs
    of the previous plotting routines.

    Inputs :
                start_val : number of first plot to use
                end_val : number of last plot to use
                folder_name : name of output folder to get results from

    Outputs :
                Saves an mp4 animation to the results folder
    '''

    from matplotlib import animation
    import matplotlib.image as mgimg

    #set up the figure
    fig = plt.figure()
    ax = plt.gca()

    #initialization of animation, plot array of zeros
    def init():
        imobj.set_data(np.zeros((250, 400)))

        return  imobj,

    def animate(i):
        ## Read in picture
        fname = os.getcwd() + '/' + folder_name + '/figs/output%d.png' %i

        ## [-1::-1], to invert the array
        # Otherwise it plots up-side down
        img = mgimg.imread(fname)[-1::-1]
        imobj.set_data(img)

        return  imobj,


    ## create an AxesImage object
    imobj = ax.imshow( np.zeros((250, 400)), origin='lower', alpha=1.0,
                       zorder=1, aspect=1.5 )
    ax.tick_params(labelbottom=False,labelleft=False)
    ax.tick_params(axis='x',bottom=False)
    ax.tick_params(axis='y',left=False)

    anim = animation.FuncAnimation(fig, animate, init_func=init, repeat = True,
                                   frames=range(start_val,end_val), interval=250,
                                   blit=True, repeat_delay=1000)

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    # fps previously 10 for the longer runs
    writer = Writer(fps=5, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'],
                    metadata=dict(artist='Me'), bitrate=-1)

    anim.save(os.getcwd()+'/' + folder_name + '/animation.mp4', writer=writer, dpi=300)
