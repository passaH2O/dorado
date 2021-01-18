# -*- coding: utf-8 -*-
"""
Higher level methods for routing and visualizing tracer particles.

Project Homepage: https://github.com/passaH2O/dorado
"""
from __future__ import division, print_function, absolute_import
from builtins import range
from .particle_track import Particles
from .particle_track import modelParams
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
def steady_plots(particle, num_iter,
                 folder_name=None, save_output=True):
    """Automated particle movement in steady flow field.

    Function to automate plotting of particle movement over a steady flow
    fields. Particles all have same number of iterations and are allowed to
    have different travel times.

    **Inputs** :

        particle : :obj:`dorado.particle_track.Particles`
            An initialized :obj:`particle_track.Particles` object with some
            generated particles.

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

    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration()
        if save_output:
            x0, y0, t0 = get_state(walk_data, 0)
            xi, yi, ti = get_state(walk_data)

            fig = plt.figure(dpi=200)
            ax = fig.add_subplot(111)
            im = ax.imshow(particle.depth)
            plt.title('Depth - Particle Iteration ' + str(i))
            cax = fig.add_axes([ax.get_position().x1+0.01,
                                ax.get_position().y0,
                                0.02,
                                ax.get_position().height])
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('Water Depth [m]')
            ax.scatter(y0, x0, c='b', s=0.75)
            ax.scatter(yi, xi, c='r', s=0.75)
            plt.savefig(folder_name+os.sep +
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


def unsteady_plots(dx, Np_tracer, seed_xloc, seed_yloc, num_steps, timestep,
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

        dx : `float`
            Length of a cell face

        Np_tracer : `int`
            Number of particles to generate.

        seed_xloc : `list`
            List of x-coordinates over which to initially distribute the
            particles.

        seed_yloc : `list`
            List of y-coordinates over which to initially distribute the
            particles.

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
    # init params
    params = modelParams()
    params.dx = dx
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

        # then define the particles class and continue
        particle = Particles(params)
        # generate some particles
        if i == 0:
            particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
        else:
            particle.generate_particles(0, [], [], 'random', walk_data)

        walk_data = particle.run_iteration(target_time=target_times[i])

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


def time_plots(particle, num_iter, folder_name=None):
    """Steady flow plots with particle travel times visualized.

    Make plots with each particle's travel time visualized.
    Routine assumes a steady flow field, but could be expanded to an unsteady
    case.

    **Inputs** :

        particle : :obj:`dorado.particle_track.Particles`
            An initialized :obj:`particle_track.Particles` object with some
            generated particles.

        num_iter : `int`
            Number of iterations to move particles

        folder_name : `str`, optional
            Path to folder in which to save output plots.

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        Saves plots and data for each iteration

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

    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration()

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
        im = ax.imshow(particle.depth)
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
    iter_exceeds_warning = 0
    # Pull out the specified value
    for ii in list(range(Np_tracer)):
        try:
            xinds.append(walk_data['xinds'][ii][iteration])
            yinds.append(walk_data['yinds'][ii][iteration])
            times.append(walk_data['travel_times'][ii][iteration])
        except IndexError:
            # If target iter exceeds walk history, return last iter
            xinds.append(walk_data['xinds'][ii][-1])
            yinds.append(walk_data['yinds'][ii][-1])
            times.append(walk_data['travel_times'][ii][-1])
            iter_exceeds_warning += 1
    
    if iter_exceeds_warning > 0:
        print('Note: %s particles have not reached %s iterations' % \
              (iter_exceeds_warning, iteration))

    return xinds, yinds, times


def get_time_state(walk_data, target_time):
    """Pull walk_data values nearest to a specific time.

    Routine to return slices of the walk_data dict at a given travel time.
    This function provides a shortcut to 'smart' indexing of the dict.

    **Inputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times

        target_time : `float`
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
        times_ii = np.array(walk_data['travel_times'][ii])
        # Find iteration nearest to target_time
        tt = np.argmin(np.abs(times_ii - target_time))

        xinds.append(walk_data['xinds'][ii][tt])
        yinds.append(walk_data['yinds'][ii][tt])
        times.append(walk_data['travel_times'][ii][tt])

        if times_ii[-1] < target_time:
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


def draw_travel_path(grid, walk_data,
                     particles_to_follow='all',
                     output_file='travel_paths.png',
                     interval=1, plot_legend=False):
    """Make a plot with the travel path of specified particles drawn out.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are
            `dorado.particle_track.modelParams.depth`,
            `dorado.particle_track.modelParams.stage`,
            `dorado.particle_track.modelParams.topography`.

        walk_data : `dict`
            Output of `steady_plots`, `unsteady_plots`, `time_plots`, as well
            as the `particle_track.run_iteration` method.

        particles_to_follow : `list`, optional
            List of particle numbers for which to draw the travel paths.
            Default is `all` particles.

        output_file : `str`, optional
            Path to save the output image to.

        interval : `int`, optional
            Determines how "smooth" each particle trajectory appears by
            skipping iterations between successive points on the path.
            Default it to show every iteration

        plot_legend : `bool`, optional
            Controls whether resulting plot includes a legend of particle IDs. 
            Default is False

    **Outputs** :

        Saves a plot of particle travel paths.

    """
    import matplotlib.patheffects as path_effects
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    if(particles_to_follow == 'all'):
        Np_tracer = len(walk_data['xinds'])
        particles_to_follow = list(range(Np_tracer))

    plt.figure(figsize=(7, 4), dpi=300)
    ax = plt.gca()
    im = ax.imshow(grid, cmap='bone', alpha=0.9)
    ax.set_title('Particle Paths')
    
    paths = [] # Place to store particle paths
    colors = [] # Place to store colors
    for i in tqdm(particles_to_follow):
        # Set color for this particle
        c = np.random.rand(3,)
        colors.append([c[0], c[1], c[2], 0.9])
        
        x = walk_data['xinds'][i][0::interval]
        y = walk_data['yinds'][i][0::interval]
        lineseg = list(zip(y, x))
        paths.append(lineseg)

    # Add new line collection to our figure, apply a background shadow
    lc = LineCollection(paths, colors=colors, 
                        linewidths=1.2, capstyle='round',
                        path_effects=[path_effects.SimpleLineShadow(offset=(0.5,-0.5),
                                                                    alpha=0.2,
                                                                    linewidth=1.6),
                                      path_effects.Normal()])
    ax.add_collection(lc)
    if plot_legend:
        custom_lines = [Line2D([0], [0], color=c, lw=1.2) for c in colors]
        ax.legend(custom_lines, [str(i) for i in particles_to_follow],
                  loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small',
                  title='Particle ID')
    plt.axis('scaled')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')


def plot_state(grid, walk_data, iteration=-1, target_time=None, c='b'):
    """Plot particle positions on an array.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are
            `dorado.particle_track.modelParams.depth`,
            `dorado.particle_track.modelParams.stage`,
            `dorado.particle_track.modelParams.topography`.

        walk_data : `dict`
            The dictionary with the particle information. This is the output
            from one of the other routines or the
            :obj:`dorado.particle_track.Particles.run_iteration()` function.

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


def snake_plots(grid,
                walk_data,
                num_steps,
                folder_name=None,
                interval=4,
                tail_length=12,
                rgba_start=[1, 0.4, 0.2, 1],
                rgba_end=[1, 0.3, 0.1, 0]):
    """Plot particle positions with a trailing tail
    
    Loops through existing walk_data history and creates a series of
    plots of the particle locations, with recent locations shown as a
    trailing tail. Default is for the tails to fade away, but they 
    can also change color. Currently this function can only sync up
    particles by iteration number, not travel_times.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are
            `dorado.particle_track.modelParams.depth`,
            `dorado.particle_track.modelParams.stage`,
            `dorado.particle_track.modelParams.topography`.

        walk_data : `dict`
            Output of `steady_plots`, `unsteady_plots`, `time_plots`, as well
            as the `particle_track.run_iteration` method.

        num_steps : `int`
            Number of states in the particle history to plot

        folder_name : `str`, optional
            Path to folder in which to save output plots

        interval : `int`, optional
            Interval of iterations to skip over between each plot. Also
            determines how "smooth" each particle trajectory appears. 
            Default is to show every 4th iteration

        tail_length : `int`, optional
            Number of previous locations to show in the "tail" behind
            each particle. Oldest location will be interval*tail_length
            iterations before the current iteration. Default is 12

        rgba_start : `list`, optional
            List of floats to use as the "recent" color in the particle
            path, specified as [red, green, blue, alpha] in the range 0-1.
            Default is a light orange

        rgba_end : `list`, optional
            List of floats to use as the "oldest" color in the particle
            path, specified as [red, green, blue, alpha] in the range 0-1.
            Default slowly fades to red-ish with an alpha of zero

    **Outputs** :

        Saves a plot of the particle travel history at each step as a png

    """
    import matplotlib.patheffects as path_effects
    from matplotlib.collections import LineCollection

    # Handle directories
    if folder_name is None:
        folder_name = os.getcwd()
    if os.path.exists(folder_name):
        print('Saving files in existing directory')
    else:
        os.makedirs(folder_name)
    if not os.path.exists(folder_name+os.sep+'figs'):
        os.makedirs(folder_name+os.sep+'figs')

    # Colors of the tail. Linearly trends from rgba_start to _end
    # Color list repeats for each particle when plotted
    colors = np.zeros((tail_length, 4))
    for c in list(range(4)):
        colors[:,c] = np.linspace(rgba_start[c], rgba_end[c],
                                  tail_length).T

    # Loop through specified number of steps
    for ii in list(range(interval, num_steps*interval, interval)):
        paths = [] # Initialize place to store paths
        # Recent chunk of the path we're going to plot
        chunks = list(range(ii, ii-tail_length*interval, -1*interval))

        # Create figure for this step
        fig = plt.figure(figsize=(7, 4), dpi=300)
        ax = plt.gca()
        im = ax.imshow(grid, cmap='bone', alpha=0.9)

        # Loop through particles
        for jj in list(range(len(walk_data['xinds']))):
            # Grab this particle's walk history
            x = walk_data['xinds'][jj]
            y = walk_data['yinds'][jj]
            lineseg = list(zip(y, x))
            newest_segment = lineseg[0:2] # Use first segment as backup

            # Check that this particle had enough iterations
            if ii > len(lineseg):
                continue # If not, skip

            # Loop through history in reverse order and grab segments
            for c in chunks:
                # Check to make sure low index isn't below zero
                if c-interval-1 >= 0:
                    # Store the most recent segment with a valid index
                    newest_segment = lineseg[c-interval-1:c:interval]
                    paths.append(newest_segment)
                else:
                    # If indices are below zero, we've reached the beginning
                    # Just re-add the most recent valid segment again
                    paths.append(newest_segment)

        # Check that we're not still plotting after all particles are done
        if len(paths) < 1:
            break

        # Add new line collection to our figure, apply a background shadow
        lc = LineCollection(paths, colors=colors,
                            linewidths=1.2, capstyle='round',
                            path_effects=[path_effects.SimpleLineShadow(offset=(0.5,-0.5),
                                                                        alpha=0.1,
                                                                        linewidth=1.6),
                                          path_effects.Normal()])
        ax.add_collection(lc)
        ax.set_title('Iteration %s' % ii)
        plt.savefig(folder_name+os.sep+'figs'
                    +os.sep+'snakefig'+str(ii)+'.png',
                    bbox_inches='tight')
        plt.close()


def show_nourishment_area(visit_freq, grid=None, walk_data=None,
                          cmap='Reds', sigma=0.7, min_alpha=0,
                          show_seed=True, seed_color='dodgerblue'):
    """Plot a smoothed, normalized heatmap of particle visit frequency

    Function will plot the history of particle travel locations in walk_data
    as a heatmap overtop the specified grid, using the output of
    :obj:`dorado.particle_track.nourishment_area()`. Colors indicate number
    of instances in which that cell was occupied by a particle.

    **Inputs** :

        visit_freq : `numpy.ndarray`
            A 2-D grid of normalized particle visit frequencies, i.e. the
            output of the :obj:`dorado.particle_track.nourishment_area()`
            function

        grid : `numpy.ndarray`, optional
            An optional 2-D grid upon which the particles will be plotted.
            Examples of grids that might be nice to use are
            `dorado.particle_track.modelParams.depth`,
            `dorado.particle_track.modelParams.stage`,
            `dorado.particle_track.modelParams.topography`.

        walk_data : `dict`, optional
            The dictionary with the particle information, which is used to
            show the seed location. This is the output from one of the other
            routines or the
            :obj:`dorado.particle_track.Particles.run_iteration()` function.

        cmap : `str`, optional
            Name of Matplotlib colormap used for the foreground (heatmap).
            Default is 'Reds'

        sigma : `float`, optional
            Degree of spatial smoothing of the heatmap used in the
            dorado.particle_track.nourishment_area() function, only used to
            adjust the plot asthetics for low sigma's

        min_alpha : `float`, optional
            Minimum alpha value of the heatmap, representing the least
            frequented cell locations, default is full transparency

        show_seed : `bool`, optional
            Determines whether resulting plot shows a marker indicating the
            (first) seed location. Uses indices of walk_data[:][0][0]. Default
            is True, but only if 'walk_data' is also provided.

        seed_color : `str`, optional
            Name of matplotlib color used for the marker seed, if shown. Default
            is a light blue to contrast with the 'Reds' heatmap.

    **Outputs** :

        ax : `matplotlib.axes`
            A `matplotlib.axes` upon which the nourishment area is drawn

    """
    from matplotlib.colors import Normalize

    # Plot heatmap with alpha based on visit_freq
    if sigma >= 0.125: # This is just a visual trial-and-error thing
        amax = np.nanpercentile(visit_freq, 60)
    else:
        amax = np.nanpercentile(visit_freq, 30)
    alphas = Normalize(0, amax, clip=True)(visit_freq) # Normalize alphas
    alphas = np.clip(alphas, min_alpha, 1)
    colors = Normalize(np.nanmin(visit_freq), 1)(visit_freq) # Normalize colors
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(colors)
    colors[..., -1] = alphas

    # Plot figure
    if len(plt.get_fignums()) < 1:
        # Create new figure axes if none are open
        fig, ax = plt.subplots(1, 1, figsize=(5,5), dpi=300)
    else:
        ax = plt.gca() # Otherwise grab existing
    ax.set_facecolor('k') # Set facecolor black
    if grid is not None:
        # Grid background intentionally dark:
        im = ax.imshow(grid, cmap='gist_gray', vmax=np.max(grid)*3)
    # Show nourishment area
    nr = ax.imshow(colors)
    plt.title('Nourishment Area')
    if (show_seed) & (walk_data is not None):
        ax.scatter(walk_data['yinds'][0][0], walk_data['xinds'][0][0],
                   c=seed_color, edgecolors='black', s=10, linewidths=0.5)
    
    return ax


def show_nourishment_time(mean_times, grid=None, walk_data=None,
                          cmap='magma', show_colorbar=True, min_alpha=0.3,
                          show_seed=True, seed_color='dodgerblue'):
    """Plot a smoothed heatmap of mean particle visit time

    Function will plot the history of mean particle travel times in walk_data
    as a heatmap overtop the specified grid, using the output of
    :obj:`dorado.particle_track.nourishment_time()`. Colors indicate the mean
    length of time particles spent in each cell (potentially after smoothing)

    **Inputs** :

        mean_times : `numpy.ndarray`
            Array of mean occupation times in each cell, i.e. the output of
            the :obj:`dorado.particle_track.nourishment_time()` function

        grid : `numpy.ndarray`, optional
            An optional 2-D grid upon which the particles will be plotted.
            Recommended to use `dorado.particle_track.modelParams.topography`

        walk_data : `dict`, optional
            The dictionary with the particle information, which is used to
            show the seed location. This is the output from one of the other
            routines or the
            :obj:`dorado.particle_track.Particles.run_iteration()` function.

        cmap : `str`, optional
            Name of Matplotlib colormap used for the foreground (heatmap).
            Default is 'magma'

        show_colorbar : `bool`, optional
            Controls whether to plot a colorbar for mean_times, default is False

        min_alpha : `float`, optional
            Minimum alpha value of the heatmap, representing the cells which
            spent the least amount of time occupied, default is 0.3

        show_seed : `bool`, optional
            Determines whether resulting plot shows a marker indicating the
            (first) seed location. Uses indices of walk_data[:][0][0]. Default
            is True, but only if 'walk_data' is also provided.

        seed_color : `str`, optional
            Name of matplotlib color used for the marker seed, if shown. Default
            is a light blue.

    **Outputs** :

        ax : `matplotlib.axes`
            A `matplotlib.axes` upon which the nourishment times are drawn

    """
    from matplotlib.colors import Normalize
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Plot heatmap with alpha based on times
    amax = np.nanpercentile(mean_times, 20)
    alphas = Normalize(0, amax, clip=True)(mean_times) # Normalize alphas
    alphas = np.clip(alphas, min_alpha, 1)
    colors = Normalize(np.nanmin(mean_times),
                       np.nanmax(mean_times))(mean_times) # Normalize colors
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(colors)
    colors[..., -1] = alphas

    # Plot figure
    if len(plt.get_fignums()) < 1:
        # Create new figure axes if none are open
        fig, ax = plt.subplots(1, 1, figsize=(5,5), dpi=300)
    else:
        ax = plt.gca() # Otherwise grab existing
    ax.set_facecolor('k') # Set facecolor black
    if grid is not None:
        # Grid background intentionally dark:
        im = ax.imshow(grid, cmap='gist_gray', vmax=np.max(grid)*3)
    # Show nourishment times
    nt = ax.imshow(colors)
    plt.title('Nourishment Times')

    # Optionally plot colorbar
    if show_colorbar:
        # Requires a few extra steps due to how we made the heatmap
        norm = matplotlib.colors.Normalize(vmin=np.nanmin(mean_times),
                                           vmax=np.nanmax(mean_times))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        divider = make_axes_locatable(ax) # Make the right size
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(sm, cax=cax)
        cbar.set_label('Mean Nourishment Time [s]')

    # Optionally show seed location
    if (show_seed) & (walk_data is not None):
        ax.scatter(walk_data['yinds'][0][0], walk_data['xinds'][0][0],
                   c=seed_color, edgecolors='black', s=10, linewidths=0.5)

    return ax
