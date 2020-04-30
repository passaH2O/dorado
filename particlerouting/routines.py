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
    Function to automate plotting of particle movement over a steady flow
    fields. Particles all have same number of iterations and are allowed to 
	have different travel times.

    Inputs :
                params : class of parameter values for the particles
                num_iter : number of iterations to move particles over
                folder_name : string of folder name to put outputs in

    Outputs :
                Script saves result of each iteration to a folder
                with the figure for each iteration as a png and the data
                with the particle locations and travel times
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

    all_walk_data = None # Initialize list for function call
    # Iterate and save results
    for i in range(0,num_iter):
        # Do particle iterations
        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data)

        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,params.Np_tracer):
            plt.scatter(all_walk_data[1][k][0],all_walk_data[0][k][0],c='b',s=0.75)
            plt.scatter(all_walk_data[1][k][-1],all_walk_data[0][k][-1],c='r',s=0.75)
        plt.imshow(params.depth)
        plt.title('Depth - Particle Iteration ' + str(i))
        cbar = plt.colorbar()
        cbar.set_label('Water Depth [m]')
        plt.axis('scaled')
        plt.savefig(os.getcwd()+'/'+folder_name+'/figs/output'+str(i)+'.png')
        plt.close()

    # save data
    np.savez(os.getcwd() + '/' + folder_name + '/data/data.npz',
             all_walk_data = all_walk_data)


def unsteady_plots(params, num_steps, timestep,
                       output_base, output_type,
                       folder_name):
    ''' 
    Function to automate plotting of particle movement in an unsteady flow
    field (time-varying). Particles all have the same travel time at the end 
    of each timestep and are allowed to have a different number of iterations.
    Flow field variables (qx, qy, depth) are updated after each timestep. 
    Because this function makes very specific assumptions about your model 
    output files, feel free to use this it as a template and change 
    the section that updates the flow field. 

    Inputs :
                params : class of particle parameter values
                num_steps : number of model timesteps being covered
                timestep : model timestep duration (seconds)
                output_base : filepath string locating hydrodynamic output files
                output_type : filetype string of the output files. Currently 
				              accepts 'csv', 'npy', and 'npz'. Assumes filenames
				              begin with either 'depth', 'stage', 'qx', 'qy', or 
				              'data', followed by timestep information
				              (limited built-in support, may require modification)
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

    # Create lists of depth, qx, qy files in the specified output_base folder
    depthlist = [x for x in os.listdir(output_base) if x.startswith('depth')]
    stagelist = [x for x in os.listdir(output_base) if x.startswith('stage')]
    qxlist = [x for x in os.listdir(output_base) if x.startswith('qx')]
    qylist = [x for x in os.listdir(output_base) if x.startswith('qy')]
    datalist = [x for x in os.listdir(output_base) if x.startswith('data')]
    if num_steps > max(len(depthlist),len(datalist)):
        print('Warning: num_steps exceeds number of model outputs in output_base')
        print('Setting num_steps to equal number of model outputs')
        num_steps = max(len(depthlist),len(datalist))

    # Create vector of target times
    target_times = np.arange(timestep, timestep*(num_steps + 1), timestep)
    all_walk_data = None
    # Iterate through model timesteps
    for i in range(0, num_steps+1):
        # load depth, stage, qx, qy for this timestep
        # Making assumption that other variables are constant between output files !!!!
        if output_type == 'csv':
            params.depth = np.loadtxt(depthlist[i], delimiter=',')
            params.stage = np.loadtxt(stagelist[i], delimiter=',')
            params.qx = np.loadtxt(qxlist[i], delimiter=',')
            params.qy = np.loadtxt(qylist[i], delimiter=',')
        elif output_type == 'npy':
            params.depth = np.load(depthlist[i])
            params.stage = np.loadtxt(stagelist[i])
            params.qx = np.load(qxlist[i])
            params.qy = np.load(qylist[i])
        elif output_type == 'npz':
            data = np.load(datalist[i])
            params.depth = data['depth']
            params.stage = data['stage']
            params.qx = data['qx']
            params.qy = data['qy']
        else:
            raise ValueError('Output datatype/structure unsupported, modify the output reading portion of the code')

        # then define the particle class and continue
        particle = Particle(params)

        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data, target_time=target_times[i])

        # make and save plots and data
        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,params.Np_tracer):
            plt.scatter(all_walk_data[1][k][0],all_walk_data[0][k][0],c='b',s=0.75)
            plt.scatter(all_walk_data[1][k][-1],all_walk_data[0][k][-1],c='r',s=0.75)
        plt.imshow(params.depth)
        plt.axis('scaled')
        plt.title('Depth at Time ' + str(target_times[i]))
        cbar = plt.colorbar()
        cbar.set_label('Water Depth [m]')
        plt.savefig(os.getcwd() + '/' + folder_name + '/figs/output'+str(i)+'.png')
        plt.close()

    # save data
    np.savez(os.getcwd() + '/' + folder_name + '/data/data.npz',
             all_walk_data = all_walk_data)


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
		
    all_walk_data = None # Initialize list for function call
    # Iterate and save results
    for i in range(0,num_iter):
        # Do particle iterations
        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data)
	
        cm = matplotlib.cm.colors.Normalize(vmax=np.max(all_walk_data[2][0:][-1]), 
                                            vmin=np.min(all_walk_data[2][0:][-1]))
        plt.figure(figsize=(4,4),dpi=200)
        for k in range(0,params.Np_tracer):
            plt.scatter(all_walk_data[1][k][0], all_walk_data[0][k][0], c='b', s=0.75)
            plt.scatter(all_walk_data[1][k][-1], all_walk_data[0][k][-1],
                        c=all_walk_data[2][k][-1],s=0.75, cmap='coolwarm', norm=cm)
        cbar = plt.colorbar()
        cbar.set_label('Particle Travel Times [s]')
        plt.imshow(params.depth)
        plt.title('Depth - Particle Iteration ' + str(i))
        cbar2 = plt.colorbar(orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.axis('scaled')
        plt.savefig(os.getcwd()+'/'+folder_name+'/figs/output'+str(i)+'.png')
        plt.close()

    # save data
    np.savez(os.getcwd() + '/' + folder_name + '/data/data.npz',
             all_walk_data = all_walk_data)


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
