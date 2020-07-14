# -*- coding: utf-8 -*-
"""
Scripts to simplify process of routing and visualizing tracer particles.

Project Homepage: https://github.com/
"""
from __future__ import division, print_function, absolute_import
from builtins import range, map
from .particle_track import Particle
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
import sys
import os
import re
import string
from tqdm import tqdm
import json


def steady_plots(params, num_iter, folder_name, save_output=True):
    """Automated particle movement in steady flow field.

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

        save_output : `bool`, optional
            Controls whether or not the output images/data are saved to disk.
            Default value is True.

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times, with
            details same as input previous_walk_data

        If `save_output` is set to True, script saves result of each iteration
        to a folder with the figure for each iteration as a png and the data
        with the particle locations and travel times as a json
        ('human-readable') text file.

    """
    # define the particle
    particle = Particle(params)

    # make directory to save the data
    if save_output:
        try:
            os.makedirs(os.getcwd() + '/' + folder_name)
            os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
            os.makedirs(os.getcwd() + '/' + folder_name + '/data')
        except Exception:
            print('Directories already exist')

    walk_data = None  # Initialize object for function call

    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration(previous_walk_data=walk_data)
        if save_output:
            fig = plt.figure(dpi=200)
            ax = plt.gca()
            im = ax.imshow(params.depth)
            plt.title('Depth - Particle Iteration ' + str(i))
            cax = fig.add_axes([ax.get_position().x1+0.01,
                                ax.get_position().y0,
                                0.02,
                                ax.get_position().height])
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('Water Depth [m]')

            for k in list(range(0, params.Np_tracer)):
                ax.scatter(walk_data['yinds'][k][0],
                           walk_data['xinds'][k][0],
                           c='b',
                           s=0.75)
                ax.scatter(walk_data['yinds'][k][-1],
                           walk_data['xinds'][k][-1],
                           c='r',
                           s=0.75)

            plt.savefig(os.getcwd()+'/'+folder_name +
                        '/figs/output'+str(i)+'.png')
            plt.close()

    if save_output:
        # save data as json text file - technically human readable
        fpath = os.getcwd() + '/' + folder_name + '/data/data.txt'
        json.dump(walk_data, open(fpath, 'w'))

        # example code to load the dictionary back from the output json file
        # data = json.load(open('data.txt'))

    return walk_data


def unsteady_plots(params, num_steps, timestep,
                   output_base, output_type,
                   folder_name):
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

        folder_name : `str`
            String of the desired output folder name

    **Outputs** :

        walk_data : `dict`
            Dictionary of all x and y locations and travel times, with
            details same as input previous_walk_data

        Script saves result of each iteration to a folder with the figure
        for each iteration as a png and the data with the particle locations
        and travel times as a json text file.

    """
    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except Exception:
        print('Directories already exist')

    # Create lists of depth, qx, qy files in the specified output_base folder
    depthlist = [x for x in os.listdir(output_base) if x.startswith('depth')]
    stagelist = [x for x in os.listdir(output_base) if x.startswith('stage')]
    qxlist = [x for x in os.listdir(output_base) if x.startswith('qx')]
    qylist = [x for x in os.listdir(output_base) if x.startswith('qy')]
    datalist = [x for x in os.listdir(output_base) if x.startswith('data')]
    if num_steps > max(len(depthlist), len(datalist)):
        print('Warning: num_steps exceeds number of model outputs in'
              ' output_base')
        print('Setting num_steps to equal number of model outputs')
        num_steps = max(len(depthlist), len(datalist))

    # Create vector of target times
    target_times = np.arange(timestep, timestep*(num_steps + 1), timestep)
    walk_data = None
    # Iterate through model timesteps
    for i in tqdm(list(range(0, num_steps+1)), ascii=True):
        # load depth, stage, qx, qy for this timestep
        # Making assumption that other variables are constant between output
        # files !!!!
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
            raise ValueError('Output datatype/structure unsupported, modify'
                             ' the output reading portion of the code')

        # then define the particle class and continue
        particle = Particle(params)

        walk_data = particle.run_iteration(previous_walk_data=walk_data,
                                           target_time=target_times[i])

        # make and save plots and data
        fig = plt.figure(dpi=200)
        for k in range(0, params.Np_tracer):
            plt.scatter(walk_data['yinds'][k][0],
                        walk_data['xinds'][k][0],
                        c='b',
                        s=0.75)
            plt.scatter(walk_data['yinds'][k][-1],
                        walk_data['xinds'][k][-1],
                        c='r',
                        s=0.75)
        ax = plt.gca()
        im = ax.imshow(params.depth)
        plt.title('Depth at Time ' + str(target_times[i]))
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,
                            0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Water Depth [m]')
        plt.savefig(os.getcwd()+'/'+folder_name+'/figs/output'+str(i)+'.png')
        plt.close()

    # save data as a json text file - technically human readable
    fpath = os.getcwd() + '/' + folder_name + '/data/data.txt'
    json.dump(walk_data, open(fpath, 'w'))

    # example code for loading the dict back from the output json file:
    # data = json.load(open('data.txt'))

    return walk_data


def time_plots(params, num_iter, folder_name):
    """Steady flow plots with particle travel times visualized.

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

        walk_data : `dict`
            Dictionary of all x and y locations and travel times, with
            details same as input previous_walk_data

        Saves plots and data for each iteration

    """
    # define the particle
    particle = Particle(params)

    # make directory to save the data
    try:
        os.makedirs(os.getcwd() + '/' + folder_name)
        os.makedirs(os.getcwd() + '/' + folder_name + '/figs')
        os.makedirs(os.getcwd() + '/' + folder_name + '/data')
    except Exception:
        print('Directories already exist')

    walk_data = None  # Initialize list for function call
    # Iterate and save results
    for i in tqdm(list(range(0, num_iter)), ascii=True):
        # Do particle iterations
        walk_data = particle.run_iteration(previous_walk_data=walk_data)

        # collect latest travel times
        temptimes = []
        for ii in list(range(0, particle.Np_tracer)):
            temptimes.append(walk_data['travel_times'][ii][-1])

        # set colorbar using 10th and 90th percentile values
        cm = matplotlib.cm.colors.Normalize(vmax=np.percentile(temptimes, 90),
                                            vmin=np.percentile(temptimes, 10))

        fig = plt.figure(dpi=200)
        plt.title('Depth - Particle Iteration ' + str(i))
        for k in range(0, params.Np_tracer):
            plt.scatter(walk_data['yinds'][k][0],
                        walk_data['xinds'][k][0],
                        c='b',
                        s=0.75)
            plt.scatter(walk_data['yinds'][k][-1],
                        walk_data['xinds'][k][-1],
                        c=walk_data['travel_times'][k][-1],
                        s=0.75,
                        cmap='coolwarm',
                        norm=cm)
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(cax=cax)
        cbar.set_label('Particle Travel Times [s]')
        im = ax.imshow(params.depth)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        cbar2 = plt.colorbar(im, cax=cax, orientation='horizontal')
        cbar2.set_label('Water Depth [m]')
        plt.savefig(os.getcwd()+'/'+folder_name+'/figs/output'+str(i)+'.png')
        plt.close()

    # save data as a json text file - technically human readable
    fpath = os.getcwd() + '/' + folder_name + '/data/data.txt'
    json.dump(walk_data, open(fpath, 'w'))

    # example code for loading the dict back from the output json file:
    # data = json.load(open('data.txt'))

    return walk_data


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
    ax = plt.gca()

    # initialization of animation, plot array of zeros
    def init():
        imobj.set_data(np.zeros((250, 400)))

        return imobj,

    def animate(i):
        # Read in picture
        fname = os.getcwd() + '/' + folder_name + '/figs/output%d.png' % i

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

    anim.save(os.getcwd()+'/' + folder_name + '/animation.mp4', writer=writer,
              dpi=300)


def draw_travel_path(depth, walk_data,
                     particles_to_follow, output_file='travel_paths.png'):
    """Make a plot with the travel path of specified particles drawn out.

    **Inputs** :

        depth : `numpy.ndarray`
            Water depth array.

        walk_data : `dict`
            Output of `steady_plots`, `unsteady_plots`, `time_plots`, as well
            as the `particle_track.run_iteration` method.

        particles_to_follow : `list`
            List of particle numbers to draw the travel paths for.

        output_file : `str`
            Path to save the output image to.

    **Outputs** :

        Saves a plot of particle travel paths.

    """
    from matplotlib import cm
    color_index = 0

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
    plt.savefig(output_file)
    plt.close()


def plot_initial(grid, all_walk_data, c='b'):
    """Plot initial particle positions on an array.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are `params.depth`, `params.stage`,
            `params.topography`.

        all_walk_data : `dict`
            The dictionary with the particle information. This is the output
            from one of the other routines or the
            :obj:particle_track.run_iteration() function.

        c : `str`, optional
            String to specify the color of the particle marks to draw on the
            figure, default is 'b' for blue

    **Outputs** :

        ax : `matplotlib.axes`
            A `matplotlib.axes` with the intended plot drawn on it

    """
    ax = plt.gca()
    ax.imshow(grid)
    x = []
    y = []
    for i in range(0, len(all_walk_data['xinds'])):
        x.append(all_walk_data['yinds'][i][0])
        y.append(all_walk_data['xinds'][i][0])
    # plot them up
    plt.scatter(x, y, c=c)

    return ax


def plot_final(grid, all_walk_data, c='r'):
    """Plot final particle positions on an array.

    **Inputs** :

        grid : `numpy.ndarray`
            A 2-D grid upon which the particles will be plotted. Examples of
            grids that might be nice to use are `params.depth`, `params.stage`,
            `params.topography`.

        all_walk_data : `dict`
            The dictionary with the particle information. This is the output
            from one of the other routines or the
            :obj:particle_track.run_iteration() function.

        c : `str`, optional
            String to specify the color of the particle marks to draw on the
            figure, default is 'r' for blue

    **Outputs** :

        ax : `matplotlib.axes`
            A `matplotlib.axes` with the intended plot drawn on it

    """
    ax = plt.gca()
    ax.imshow(grid)
    x = []
    y = []
    for i in range(0, len(all_walk_data['xinds'])):
        x.append(all_walk_data['yinds'][i][-1])
        y.append(all_walk_data['xinds'][i][-1])

    plt.scatter(x, y, c=c)

    return ax
