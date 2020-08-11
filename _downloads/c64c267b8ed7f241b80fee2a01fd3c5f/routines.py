# -*- coding: utf-8 -*-
"""
Higher level methods for routing and visualizing tracer particles.

Project Homepage: https://github.com/passaH2O/dorado
"""
from __future__ import division, print_function, absolute_import
from builtins import range
from .particle_track import Particle
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
import os
from tqdm import tqdm
import json


# --------------------------------------------------------
# Functions for running random walk model
# --------------------------------------------------------
def steady_plots(params, num_iter, folder_name=None, save_output=True):
    """Automated particle movement in steady flow field.

    Function to automate plotting of particle movement over a steady flow
    fields. Particles all have same number of iterations and are allowed to
    have different travel times.

    **Inputs** :

        params : `obj`
            Class of parameter values for the particles

        num_iter : `int`
            Number of iterations to move particles over

        folder_name : `str`, optional
            Path to folder in which to save output plots.

        save_output : `bool`, optional
            Controls whether or not the output images/data are saved to disk.
            Default value is True.

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        If `save_output` is set to True, script saves result of each iteration
        to a folder with the figure for each iteration as a png and the data
        with the particle locations and travel times as a json
        ('human-readable') text file.

    """
    # define the particle
    particle = Particle(params)

    # make directory to save the data
    if save_output:
        if folder_name is None:
            folder_name = os.getcwd()
        if os.path.exists(folder_name):
            print('Saving files in existing directory')
        else:
            os.makedirs(folder_name)
        if not os.path.exists(folder_name+os.sep+'figs'):
            os.makedirs(folder_name+os.sep+'figs')
        if not os.path.exists(folder_name+os.sep+'data'):
            os.makedirs(folder_name+os.sep+'data')

    walk_data = None  # Initialize object for function call

    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration(previous_walk_data=walk_data)
        if save_output:
            x0, y0, t0 = get_state(walk_data, 0)
            xi, yi, ti = get_state(walk_data)

            fig = plt.figure(dpi=200)
            ax = fig.add_subplot(111)
            im = ax.imshow(params.depth)
            plt.title('Depth - Particle Iteration ' + str(i))
            cax = fig.add_axes([ax.get_position().x1+0.01,
                                ax.get_position().y0,
                                0.02,
                                ax.get_position().height])
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('Water Depth [m]')
            ax.scatter(y0, x0, c='b', s=0.75)
            ax.scatter(yi, xi, c='r', s=0.75)
            plt.savefig(folder_name+os.sep+
                        'figs'+os.sep+'output'+str(i)+'.png',
                        bbox_inches='tight')
            plt.close()

    if save_output:
        # save data as json text file - technically human readable
        fpath = folder_name+os.sep+'data'+os.sep+'data.txt'
        json.dump(walk_data, open(fpath, 'w'))

        # example code to load the dictionary back from the output json file
        # data = json.load(open('data.txt'))

    return walk_data


def unsteady_plots(params, num_steps, timestep,
                   output_base, output_type,
                   folder_name=None):
    """Automated particle movement in unsteady flow.

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

        folder_name : `str`, optional
            Path to folder in which to save output plots.

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        Script saves result of each iteration to a folder with the figure
        for each iteration as a png and the data with the particle locations
        and travel times as a json text file.

    """
    # make directory to save the data
    if folder_name is None:
        folder_name = os.getcwd()
    if os.path.exists(folder_name):
        print('Saving files in existing directory')
    else:
        os.makedirs(folder_name)
    if not os.path.exists(folder_name+os.sep+'figs'):
        os.makedirs(folder_name+os.sep+'figs')
    if not os.path.exists(folder_name+os.sep+'data'):
        os.makedirs(folder_name+os.sep+'data')

    # Create lists of depth, qx, qy files in the specified output_base folder
    depthlist = [x for x in os.listdir(output_base) if x.startswith('depth')]
    stagelist = [x for x in os.listdir(output_base) if x.startswith('stage')]
    qxlist = [x for x in os.listdir(output_base) if x.startswith('qx')]
    qylist = [x for x in os.listdir(output_base) if x.startswith('qy')]
    datalist = [x for x in os.listdir(output_base) if x.startswith('data')]
    # sort the lists
    depthlist = sorted(depthlist)
    stagelist = sorted(stagelist)
    qxlist = sorted(qxlist)
    qylist = sorted(qylist)
    datalist = sorted(datalist)
    if num_steps > max(len(depthlist), len(datalist)):
        print('Warning: num_steps exceeds number of model outputs in'
              ' output_base')
        print('Setting num_steps to equal number of model outputs')
        num_steps = max(len(depthlist), len(datalist))

    # Create vector of target times
    target_times = np.arange(timestep, timestep*(num_steps + 1), timestep)
    walk_data = None
    # Iterate through model timesteps
    for i in tqdm(list(range(0, num_steps)), ascii=True):
        # load depth, stage, qx, qy for this timestep
        # Making assumption that other variables are constant between output
        # files !!!!
        if output_type == 'csv':
            params.depth = np.loadtxt(os.path.join(output_base, depthlist[i]),
                                      delimiter=',')
            params.stage = np.loadtxt(os.path.join(output_base, stagelist[i]),
                                      delimiter=',')
            params.qx = np.loadtxt(os.path.join(output_base, qxlist[i]),
                                   delimiter=',')
            params.qy = np.loadtxt(os.path.join(output_base, qylist[i]),
                                   delimiter=',')
        elif output_type == 'npy':
            params.depth = np.load(os.path.join(output_base, depthlist[i]))
            params.stage = np.load(os.path.join(output_base, stagelist[i]))
            params.qx = np.load(os.path.join(output_base, qxlist[i]))
            params.qy = np.load(os.path.join(output_base, qylist[i]))
        elif output_type == 'npz':
            data = np.load(os.path.join(output_base, datalist[i]))
            params.depth = data['depth']
            params.stage = data['stage']
            params.qx = data['qx']
            params.qy = data['qy']
        else:
            raise ValueError('Output datatype/structure unsupported, modify'
                             ' the output reading portion of the code')

        # then define the particle class and continue
        particle = Particle(params)

        walk_data = particle.run_iteration(previous_walk_data=walk_data,
                                           target_time=target_times[i])

        x0, y0, t0 = get_state(walk_data, 0)
        xi, yi, ti = get_state(walk_data)

        # make and save plots and data
        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(111)
        ax.scatter(y0, x0, c='b', s=0.75)
        ax.scatter(yi, xi, c='r', s=0.75)
        im = ax.imshow(params.depth)
        plt.title('Depth at Time ' + str(target_times[i]))
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,
                            0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Water Depth [m]')
        plt.savefig(folder_name+os.sep+'figs'+os.sep+'output'+str(i)+'.png',
                    bbox_inches='tight')
        plt.close()

    # save data as a json text file - technically human readable
    fpath = folder_name+os.sep+'data'+os.sep+'data.txt'
    json.dump(walk_data, open(fpath, 'w'))

    # example code for loading the dict back from the output json file:
    # data = json.load(open('data.txt'))

    return walk_data


def time_plots(params, num_iter, folder_name=None):
    """Steady flow plots with particle travel times visualized.

    Make plots with each particle's travel time visualized.
    Routine assumes a steady flow field, but could be expanded to an unsteady
    case.

    **Inputs** :

        params : `obj`
            Parameters for the particle

        num_iter : `int`
            Number of iterations to move particles

        folder_name : `str`, optional
            Path to folder in which to save output plots.

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        Saves plots and data for each iteration

    """
    # define the particle
    particle = Particle(params)

    # make directory to save the data
    if folder_name is None:
        folder_name = os.getcwd()
    if os.path.exists(folder_name):
        print('Saving files in existing directory')
    else:
        os.makedirs(folder_name)
    if not os.path.exists(folder_name+os.sep+'figs'):
        os.makedirs(folder_name+os.sep+'figs')
    if not os.path.exists(folder_name+os.sep+'data'):
        os.makedirs(folder_name+os.sep+'data')

    walk_data = None  # Initialize list for function call
    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration(previous_walk_data=walk_data)

        # collect latest travel times
        x0, y0, t0 = get_state(walk_data, 0)
        xi, yi, temptimes = get_state(walk_data)

        # set colorbar using 10th and 90th percentile values
        cm = matplotlib.cm.colors.Normalize(vmax=np.percentile(temptimes, 90),
                                            vmin=np.percentile(temptimes, 10))

        fig = plt.figure(dpi=200)
        ax = plt.gca()
        plt.title('Depth - Particle Iteration ' + str(i))
        ax.scatter(y0, x0, c='b', s=0.75)
        sc = ax.scatter(yi, xi, c=temptimes, s=0.75, cmap='coolwarm', norm=cm)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(sc, cax=cax)
        cbar.set_label('Particle Travel Times [s]')
        im = ax.imshow(params.depth)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        cbar2 = plt.colorbar(im, cax=cax, orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.savefig(folder_name+os.sep+'figs'+os.sep+'output'+str(i)+'.png',
                    bbox_inches='tight')
        plt.close()

    # save data as a json text file - technically human readable
    fpath = folder_name+os.sep+'data'+os.sep+'data.txt'
    json.dump(walk_data, open(fpath, 'w'))

    # example code for loading the dict back from the output json file:
    # data = json.load(open('data.txt'))

    return walk_data


# --------------------------------------------------------
# Functions for plotting/interpreting outputs
# --------------------------------------------------------
def get_state(walk_data, iteration=-1):
    """Pull walk_data values from a specific iteration.

    Routine to return slices of the walk_data dict at a given iteration #.
    This function provides a shortcut to 'smart' indexing of the dict.
    By default returns information on the most recent step

    **Inputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        iteration : `int`, optional
            Iteration number at which to slice the dictionary. Defaults
            to -1, i.e. the most recent step

    **Outputs** :

        xinds : `list`
            List containing equivalent of walk_data['xinds'][:][iteration]

        yinds : `list`
            List containing equivalent of walk_data['yinds'][:][iteration]

        times : `list`
            List containing equivalent of
            walk_data['travel_times'][:][iteration]

    """
    iteration = int(iteration)
    Np_tracer = len(walk_data['xinds'])  # Number of particles

    xinds = []
    yinds = []
    times = []
    # Pull out the specified value
    for ii in list(range(Np_tracer)):
        xinds.append(walk_data['xinds'][ii][iteration])
        yinds.append(walk_data['yinds'][ii][iteration])
        times.append(walk_data['travel_times'][ii][iteration])

    return xinds, yinds, times


def get_time_state(walk_data, target_time):
    """Pull walk_data values nearest to a specific time.

    Routine to return slices of the walk_data dict at a given travel time.
    This function provides a shortcut to 'smart' indexing of the dict.

    **Inputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        target_time : `float`, optional
            Travel time at which to slice the dictionary.

    **Outputs** :

        xinds : `list`
            List containing equivalent of walk_data['xinds'][:][i], where i
            represents the index at which travel time is nearest to input.

        yinds : `list`
            List containing equivalent of walk_data['yinds'][:][i], where i
            represents the index at which travel time is nearest to input.

        times : `list`
            List containing equivalent of walk_data['travel_times'][:][i],
            where i represents the index at which travel time is nearest
            to input. Times will differ slightly from input time due to
            the nature of the method.

    """
    Np_tracer = len(walk_data['xinds'])  # Number of particles

    xinds = []
    yinds = []
    times = []
    # Pull out the specified value
    for ii in list(range(Np_tracer)):
        for jj in list(range(len(walk_data['travel_times'][ii])-1)):
            this_time = walk_data['travel_times'][ii][jj]
            next_time = walk_data['travel_times'][ii][jj+1]
            # Save once the next step takes us farther from the target time
            # than we are right now
            if abs(this_time - target_time) <= abs(next_time - target_time):
                xinds.append(walk_data['xinds'][ii][jj])
                yinds.append(walk_data['yinds'][ii][jj])
                times.append(walk_data['travel_times'][ii][jj])
                break
            # If we made it to the end without getting there, save last
            elif (jj == len(walk_data['travel_times'][ii])-1):
                xinds.append(walk_data['xinds'][ii][jj+1])
                yinds.append(walk_data['yinds'][ii][jj+1])
                times.append(walk_data['travel_times'][ii][jj+1])
                print('Note: Particle '+str(ii)+' never reached target_time')

    return xinds, yinds, times


def plot_exposure_time(walk_data,
                       exposure_times,
                       folder_name=None,
                       timedelta=1,
                       nbins=100,
                       save_output=True):
    """Plot exposure time distribution of particles in a specified region.

    Function to plot the exposure time distribution (ETD) of particles to
    the specified region. Relies on the output of particle_track.exposure_time

    **Inputs** :

        walk_data : `dict`
            Output of a previous function call to run_iteration.

        exposure_times : `list`
            List [] of floats containing the output of
            particle_track.exposure_time

        folder_name : `str`, optional
            Path to folder in which to save output plots.

        timedelta : `int or float`, optional
            Unit of time for time-axis of ETD plots, specified as time
            in seconds (e.g. an input of 60 plots things by minute).

        nbins : `int`, optional
            Number of bins to use as the time axis for differential ETD.
            Using fewer bins smoothes out curves.

        save_output : `bool`, optional
            Controls whether or not the output images/data are saved to disk.
            Default value is True.

    **Outputs** :

        If `save_output` is set to True, script saves plots of the cumulative
        and differential forms of the ETD as a png and the list of exposure
        times as a json ('human-readable') text file.
    """
    # Initialize arrays to record exposure time of each particle
    Np_tracer = len(walk_data['xinds'])  # Number of particles
    # Record final travel times
    x, y, end_time = get_state(walk_data)

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

    # Handle directories
    if folder_name is None:
        folder_name = os.getcwd()
    if os.path.exists(folder_name):
        print('Saving files in existing directory')
    else:
        os.makedirs(folder_name)
    if not os.path.exists(folder_name+os.sep+'figs'):
        os.makedirs(folder_name+os.sep+'figs')
    if not os.path.exists(folder_name+os.sep+'data'):
        os.makedirs(folder_name+os.sep+'data')

    # Save exposure times by particle ID
    fpath = folder_name+os.sep+'data'+os.sep+'exposure_times.txt'
    json.dump(exposure_times, open(fpath, 'w'))
    exposure_times = np.array(exposure_times)

    # Set end of ETD as the minimum travel time of particles
    # Exposure times after that are unreliable because not all particles have
    # traveled for that long
    end_time = min(end_time)

    # Ignore particles that never entered ROI or exited ROI for plotting
    # If never entered, ET of 0
    plotting_times = exposure_times[exposure_times > 1e-6]
    # If never exited, ET ~= end_time
    plotting_times = plotting_times[plotting_times < 0.99*end_time]
    # Number of particles included in ETD plots
    num_particles_included = len(plotting_times)

    # Full time vector (x-values) of CDF
    # Add origin for plot
    full_time_vect = np.append([0], np.sort(plotting_times))
    full_time_vect = np.append(full_time_vect, [end_time])
    # Y-values of CDF, normalized
    frac_exited = np.arange(0, num_particles_included + 1,
                            dtype='float')/Np_tracer
    frac_exited = np.append(frac_exited,
                            [float(num_particles_included)/float(Np_tracer)])

    if save_output:
        # Plot the cumulative ETD in its exact form
        plt.figure(figsize=(5, 3), dpi=150)
        plt.step(full_time_vect/timedelta, frac_exited, where='post')
        plt.title('Cumulative Exposure Time Distribution')
        plt.xlabel('Time ' + timeunit)
        plt.ylabel('F(t) [-]')
        plt.xlim([0, end_time/timedelta])
        plt.ylim([0, 1])
        plt.savefig(folder_name+os.sep+'figs'+os.sep+'Exact_CETD.png',
                    bbox_inches='tight')
        plt.close()

    # Smooth out the CDF by making it regular in time.
    # Here we use 'previous' interpolation to be maximally accurate in time
    create_smooth_CDF = scipy.interpolate.interp1d(full_time_vect,
                                                   frac_exited,
                                                   kind='previous')
    smooth_time_vect = np.linspace(0, end_time, nbins)
    smooth_CDF = create_smooth_CDF(smooth_time_vect)

    # Plot the cumulative ETD in its smooth form
    plt.figure(figsize=(5, 3), dpi=150)
    plt.plot(smooth_time_vect/timedelta, smooth_CDF)
    plt.title('Cumulative Exposure Time Distribution')
    plt.xlabel('Time ' + timeunit)
    plt.ylabel('F(t) [-]')
    plt.xlim([0, end_time/timedelta])
    plt.ylim([0, 1])
    if save_output:
        plt.savefig(folder_name+os.sep+'figs'+os.sep+'Smooth_CETD.png',
                    bbox_inches='tight')

    # Derive differential ETD from the CDF
    # Here we use 'linear' interpolation, because 'previous'
    # produces a choppy derivative if there aren't enough particles
    create_linear_CDF = scipy.interpolate.interp1d(full_time_vect,
                                                   frac_exited,
                                                   kind='linear')
    linear_CDF = create_linear_CDF(smooth_time_vect)

    timestep = smooth_time_vect[1] - smooth_time_vect[0]
    RTD = np.gradient(linear_CDF, timestep)

    plt.figure(figsize=(5, 4), dpi=150)
    plt.plot(smooth_time_vect/timedelta, RTD*timedelta)
    plt.title('Exposure Time Distribution')
    plt.xlabel('Time ' + timeunit)
    plt.ylabel('E(t) ' + timeunit[0:-1] + '$^{-1}$]')
    plt.xlim([0, end_time/timedelta])
    plt.ylim(ymin=0)
    if save_output:
        plt.savefig(folder_name+os.sep+'figs'+os.sep+'ETD.png',
                    bbox_inches='tight')


# Function to automate animation of the png outputs
# requires installation of the animation writer 'ffmpeg' which is not part of
# the default installation set of packages
def animate_plots(start_val, end_val, folder_name):
    """Animation routine, requires an optional animation writer.

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

    """
    from matplotlib import animation
    import matplotlib.image as mgimg

    # set up the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # initialization of animation, plot array of zeros
    def init():
        imobj.set_data(np.zeros((250, 400)))

        return imobj,

    def animate(i):
        # Read in picture
        fname = folder_name+os.sep+'figs'+os.sep+'output%d.png' % i

        # [-1::-1], to invert the array
        # Otherwise it plots up-side down
        img = mgimg.imread(fname)[-1::-1]
        imobj.set_data(img)

        return imobj,

    # create an AxesImage object
    imobj = ax.imshow(np.zeros((250, 400)), origin='lower', alpha=1.0,
                      zorder=1, aspect=1.5)
    ax.tick_params(labelbottom=False, labelleft=False)
    ax.tick_params(axis='x', bottom=False)
    ax.tick_params(axis='y', left=False)

    anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=True,
                                   frames=list(range(start_val, end_val)),
                                   interval=250,
                                   blit=True, repeat_delay=1000)

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    # fps previously 10 for the longer runs
    writer = Writer(fps=5, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'],
                    metadata=dict(artist='Me'), bitrate=-1)

    anim.save(folder_name+os.sep+'animation.mp4', writer=writer,
              dpi=300)


def draw_travel_path(depth, walk_data,
                     particles_to_follow='all',
                     output_file='travel_paths.png'):
    """Make a plot with the travel path of specified particles drawn out.

    **Inputs** :

        depth : `numpy.ndarray`
            Water depth array.

        walk_data : `dict`
            Output of `steady_plots`, `unsteady_plots`, `time_plots`, as well
            as the `particle_track.run_iteration` method.

        particles_to_follow : `list`, optional
            List of particle numbers for which to draw the travel paths.
            Default is `all` particles.

        output_file : `str`, optional
            Path to save the output image to.

    **Outputs** :

        Saves a plot of particle travel paths.

    """
    from matplotlib import cm
    color_index = 0

    if(particles_to_follow == 'all'):
        Np_tracer = len(walk_data['xinds'])
        particles_to_follow = list(range(Np_tracer))

    plt.figure(figsize=(7, 4), dpi=300)
    plt.imshow(depth, cmap='bone')
    plt.title('Particle paths overlaid on water depth')
    cbar2 = plt.colorbar(orientation='horizontal')
    cbar2.set_label('Water Depth [m]')

    for i in tqdm(particles_to_follow):
        # set color for this particle (using a discrete colormap)
        c = cm.Set1(color_index)
        col = [c[0], c[1], c[2], 0.85]  # make color a bit transparent
        color_index += 1
        # visualize this particle's travel path
        for j in range(1, len(walk_data['xinds'][i][:])):
            # define old x-y point
            old_x = walk_data['xinds'][i][j-1]
            old_y = walk_data['yinds'][i][j-1]
            # identify new x-y point
            new_x = walk_data['xinds'][i][j]
            new_y = walk_data['yinds'][i][j]
            # add the line for the particle
            if j == 1:
                plt.plot([old_y, new_y],
                         [old_x, new_x],
                         color=col,
                         linewidth=0.9,
                         label='Particle ' + str(i))

            else:
                plt.plot([old_y, new_y],
                         [old_x, new_x],
                         color=col,
                         linewidth=0.9,
                         label='_nolegend_')

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.axis('scaled')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()


def plot_state(grid, walk_data, iteration=-1, target_time=None, c='b'):
    """Plot particle positions on an array.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are `params.depth`, `params.stage`,
            `params.topography`.

        walk_data : `dict`
            The dictionary with the particle information. This is the output
            from one of the other routines or the
            :obj:particle_track.run_iteration() function.

        iteration : `int`, optional
            Iteration number at which to plot the particles. Default is -1,
            i.e. the most recent step

        target_time : `float`, optional
            Travel time at which to plot the particles. If specified, iteration
            is ignored.

        c : `str`, optional
            String to specify the color of the particle marks to draw on the
            figure, default is 'b' for blue

    **Outputs** :

        ax : `matplotlib.axes`
            A `matplotlib.axes` with the intended plot drawn on it

    """
    if target_time is None:
        x, y, t = get_state(walk_data, int(iteration))
    else:
        x, y, t = get_time_state(walk_data, target_time)
    ax = plt.gca()
    ax.imshow(grid)
    plt.scatter(y, x, c=c)

    return ax
