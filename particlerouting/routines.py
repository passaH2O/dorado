# -*- coding: utf-8 -*-
"""
Scripts to simplify process of routing and visualizing tracer particles.

Project Homepage: https://github.com/
"""
from __future__ import division, print_function, absolute_import
from builtins import range, map
from .particle_track import Particle
from math import floor, sqrt, pi
import numpy as np
from random import shuffle
import matplotlib
from matplotlib import pyplot as plt
import scipy
from scipy import ndimage
from scipy.sparse import lil_matrix, csc_matrix, hstack
import sys, os, re, string
from netCDF4 import Dataset
import time as time_lib
import logging
from tqdm import tqdm

def steady_plots(params,num_iter,folder_name):
    '''
    Function to automate plotting of particle movement over a steady flow
    fields. Particles all have same number of iterations and are allowed to
	have different travel times.

    **Inputs** :

        params : `obj`
            Class of parameter values for the particles

        num_iter : `int`
            Number of iterations to move particles over

        folder_name : `str`
            String of folder name to put outputs in

    **Outputs** :

        all_walk_data : `list`
            Nested list of all x and y locations and travel times, with
            details same as input previous_walk_data

        Script saves result of each iteration to a folder with the figure for
        each iteration as a png and the data with the particle locations and
        travel times

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
    for i in tqdm(list(range(0,num_iter)), ascii=True):
        # Do particle iterations
        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data)
        plt.figure(figsize=(4,4),dpi=200)
        for k in list(range(0,params.Np_tracer)):
            plt.scatter(all_walk_data[1][k][0],all_walk_data[0][k][0],
                        c='b',s=0.75)
            plt.scatter(all_walk_data[1][k][-1],
                        all_walk_data[0][k][-1],c='r',s=0.75)
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

    return all_walk_data


def unsteady_plots(params, num_steps, timestep,
                       output_base, output_type,
                       folder_name):
    '''
    Function to automate plotting of particle movement in an unsteady flow
    field (time-varying). Particles all have the same travel time at the end
    of each timestep and are allowed to have a different number of iterations.
    Flow field variables (qx, qy, depth) are updated after each timestep.
    Because this function makes very specific assumptions about your model
    output files, feel free to use this as a template and change
    the section that updates the flow field.

    **Inputs** :

        params : `obj`
            Class of particle parameter values

        num_steps : `int`
            Number of model timesteps being covered

        timestep : `float`
            Model timestep duration (seconds)

        output_base : `str`
            Filepath string locating hydrodynamic output files

        output_type : `str`
            Filetype string of the output files. Currently accepts 'csv',
            'npy', and 'npz'. Assumes filenames begin with either 'depth',
            'stage', 'qx', 'qy', or 'data', followed by timestep information
		    (limited built-in support, may require modification)

        folder_name : `str`
            String of the desired output folder name

    **Outputs** :

        all_walk_data : `list`
            Nested list of all x and y locations and travel times, with
            details same as input previous_walk_data

        Script saves result of each iteration to a folder with both the figure
        for each iteration as a png and the data with the particle start and
        end locations

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
    for i in tqdm(list(range(0, num_steps+1)), ascii=True):
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

    return all_walk_data


def time_plots(params,num_iter,folder_name):
    '''
    Make plots with each particle's travel time visualized.
    Routine assumes a steady flow field, but could be expanded to an unsteady
    case.

    **Inputs** :

        params : `obj`
            Parameters for the particle

        num_iter : `int`
            Number of iterations to move particles

        folder_name : `str`
            String of desired output folder name

    **Outputs** :

        all_walk_data : `list`
            Nested list of all x and y locations and travel times, with
            details same as input previous_walk_data

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
    for i in tqdm(list(range(0,num_iter)), ascii=True):
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

    return all_walk_data


### Function to automate animation of the png outputs
# requires installation of the animation writer 'ffmpeg' which is not part of the
# default PartingRouting installation set of packages
def animate_plots(start_val,end_val,folder_name):
    '''
    Routine to make mp4 animation of the particle routing from png outputs
    of the previous plotting routines.

    **Inputs** :

        start_val : `int`
            Number of first plot to use

        end_val : `int`
            Number of last plot to use

        folder_name : `str`
            Name of output folder to get results from

    **Outputs** :

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
                                   frames=list(range(start_val,end_val)), interval=250,
                                   blit=True, repeat_delay=1000)

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    # fps previously 10 for the longer runs
    writer = Writer(fps=5, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'],
                    metadata=dict(artist='Me'), bitrate=-1)

    anim.save(os.getcwd()+'/' + folder_name + '/animation.mp4', writer=writer, dpi=300)


def exposure_time(all_walk_data,
                  region_of_interest,
                  folder_name,
                  timedelta=1,
                  nbins=100):
    '''
    Routine to measure the exposure time distribution (ETD) of particles to
    the specified region. For steady flows, the ETD is exactly equivalent to
    the residence time distribution. For unsteady flows, if particles make
    multiple excursions into the region, all of those times are counted.

    **Inputs** :

        all_walk_data : `list`
            Output of a previous function call to run_iteration

        region_of_interest : `int array`
            Binary array the same size as input arrays in params class
            with 1's everywhere inside the region in which we want to
            measure exposure time, and 0's everywhere else.

        folder_name : `str`
            Name of output folder to get results from

        timedelta : `int or float`
            Unit of time for time-axis of ETD plots, specified as time
            in seconds (e.g. an input of 60 plots things by minute)

        nbins : `int`
            Number of bins to use as the time axis for differential ETD.
            Using fewer bins smooths out curves

    **Outputs** :

        exposure_times : `list`
            List of exposure times to region of interest, listed
            in order of particle ID

        Saves plots of the cumulative and differential forms of the ETD
    '''
    # Initialize arrays to record exposure time of each particle
    Np_tracer = len(all_walk_data[0]) # Number of particles
    exposure_times = np.zeros([Np_tracer], dtype='float') # Array to be populated
    end_time = np.zeros([Np_tracer], dtype='float') # Array to record final travel times

    # Handle the timedelta
    if timedelta == 1:
        timeunit = '[s]'
    elif timedelta == 60:
        timeunit = '[m]'
    elif timedelta == 3600:
        timeunit = '[hr]'
    elif timedelta == 86400:
        timeunit == '[day]'
    else:
        timeunit = '[' + str(timedelta) + ' s]'

    # Loop through particles to measure exposure time
    for ii in tqdm(list(range(0, Np_tracer)), ascii=True):
        # Determine the starting region for particle ii
        previous_reg = region_of_interest[int(all_walk_data[0][ii][0]), int(all_walk_data[1][ii][0])]
        end_time[ii] = all_walk_data[2][ii][-1] # Length of runtime for particle ii

        # Loop through iterations
        for jj in list(range(1, len(all_walk_data[2][ii]))):
            # Determine the new region and compare to previous region
            current_reg = region_of_interest[int(all_walk_data[0][ii][jj]), int(all_walk_data[1][ii][jj])]

            # Check to see if whole step was inside ROI
            if (current_reg + previous_reg) == 2: # If so, travel time of the whole step added to ET
                exposure_times[ii] += (all_walk_data[2][ii][jj] - all_walk_data[2][ii][jj-1])
            # Check to see if half of the step was inside ROI (either entering or exiting)
            elif (current_reg + previous_reg) == 1: # If so, travel time of half of the step added to ET
                exposure_times[ii] += 0.5*(all_walk_data[2][ii][jj] - all_walk_data[2][ii][jj-1])

            # Update previous region
            previous_reg = current_reg

            # Check if particle is still stuck in ROI at the end of the run, which can bias result
            if jj == len(all_walk_data[2][ii])-1:
                if current_reg == 1:
                    print(('Warning: Particle ' + str(ii) + ' is still within ROI at final timestep. \n' + \
                           'Run more iterations to get tail of ETD'))

    # Set end of ETD as the minimum travel time of particles
    # Exposure times after that are unreliable because not all particles have traveled for that long
    end_time = min(end_time)

    # Ignore particles that never entered ROI for plotting
    plotting_times = exposure_times[exposure_times > 1e-6] # Those particles will have had an ET of 0
    num_particles_included = len(plotting_times) # Number of particles that spent at least some time in ROI

    # Full time vector (x-values) of CDF
    full_time_vect = np.append([0], np.sort(plotting_times)) # Add origin for plot
    full_time_vect = np.append(full_time_vect, [end_time])
    # Y-values of CDF, normalized
    frac_exited = np.arange(0, num_particles_included + 1, dtype = 'float')/Np_tracer
    frac_exited = np.append(frac_exited, [float(num_particles_included)/float(Np_tracer)])

    # Plot the cumulative ETD in its exact form
    plt.figure(figsize=(5,3), dpi=150)
    plt.step(full_time_vect/timedelta, frac_exited, where='post')
    plt.title('Cumulative Exposure Time Distribution')
    plt.xlabel('Time ' + timeunit)
    plt.ylabel('F(t) [-]')
    plt.xlim([0, end_time/timedelta])
    plt.ylim([0, 1])
    plt.savefig(os.getcwd()+'/'+folder_name+'/figs/Exact_CETD.png')
    plt.close()

    # Smooth out the CDF by making it regular in time.
    # Here we use 'previous' interpolation to be maximally accurate in time
    create_smooth_CDF = scipy.interpolate.interp1d(full_time_vect, frac_exited, kind = 'previous')
    smooth_time_vect = np.linspace(0, end_time, nbins)
    smooth_CDF = create_smooth_CDF(smooth_time_vect)

    # Plot the cumulative ETD in its smooth form
    plt.figure(figsize=(5,3), dpi=150)
    plt.plot(smooth_time_vect/timedelta, smooth_CDF)
    plt.title('Cumulative Exposure Time Distribution')
    plt.xlabel('Time ' + timeunit)
    plt.ylabel('F(t) [-]')
    plt.xlim([0, end_time/timedelta])
    plt.ylim([0, 1])
    plt.savefig(os.getcwd()+'/'+folder_name+'/figs/Smooth_CETD.png')
    plt.close

    # Derive differential ETD from the CDF. Here we use 'linear' interpolation, because 'previous'
    # produces a choppy derivative if there aren't enough particles
    create_linear_CDF = scipy.interpolate.interp1d(full_time_vect, frac_exited, kind = 'linear')
    linear_CDF = create_linear_CDF(smooth_time_vect)

    timestep = smooth_time_vect[1] - smooth_time_vect[0]
    RTD = np.gradient(linear_CDF, timestep)

    plt.figure(figsize=(5,4), dpi=150)
    plt.plot(smooth_time_vect/timedelta, RTD*timedelta)
    plt.title('Exposure Time Distribution')
    plt.xlabel('Time ' + timeunit)
    plt.ylabel('E(t) ' + timeunit[0:-1] + '$^{-1}$]')
    plt.xlim([0, end_time/timedelta])
    plt.savefig(os.getcwd()+'/'+folder_name+'/figs/ETD.png')
    plt.close

    return exposure_times
