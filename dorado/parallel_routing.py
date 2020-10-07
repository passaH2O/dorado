"""Parallel functionality using the native Python multiprocessing library.

Parallel functionality using multiprocessing (included with Python).
For local parallelization using CPU cores.

Project Homepage: https://github.com/passaH2O/dorado
"""
from __future__ import division, print_function, absolute_import
from builtins import range, map
import numpy as np
import sys
import os
import re
import string
from multiprocessing import Pool


# convert run script into a function that returns a single dictionary 'data'
def run_iter(pobj):
    """Wrapper run script for the particle iterations.

    This is an internal function that should not be called directly. Rather it
    is wrapped by the `parallel_routing` function and mapped to the CPU cores
    to allow for particle routing in parallel. Uses a small `parallel_obj`
    class to hold a bunch of attributes and do the particle generation and
    routing.

    **Inputs** :

        pobj : :obj:`parallel_obj`
            Special parallel object.

    **Outputs** :

        all_walk_data : `list`
            Nested list of all x and y locations and travel times, with
            details same as input previous_walk_data

    """
    pobj.particles.generate_particles(pobj.Np_tracer, pobj.seed_xloc,
                                      pobj.seed_yloc)
    # do iterations
    for i in list(range(0, pobj.num_iter)):
        all_walk = pobj.particles.run_iteration()

    return all_walk


def combine_result(par_result):
    """Combine results from each core.

    Take the parallel resulting list and combine the entries so that a single
    dictionary with the beginning/end indices and travel times for all
    particles is returned as opposed to a list with one dictionary per
    parallel process

    **Inputs** :

        par_result : `list`
            List of length(num_cores) with a dictionary of the beg/end indices
            and travel times for each particle computed by that process/core

    **Outputs** :

        single_result : `list`
            Nested list that matches 'all_walk_data'

    """
    # initiate final results dictionary
    single_result = dict()
    single_result['x_inds'] = []
    single_result['y_inds'] = []
    single_result['travel_times'] = []

    # populate the dictionary
    # loop through results for each core
    for i in range(0, len(par_result)):
        # append results for each category
        for j in range(0, len(par_result[i][0])):
            single_result['x_inds'].append(par_result[i][0][j])
            single_result['y_inds'].append(par_result[i][1][j])
            single_result['travel_times'].append(par_result[i][2][j])

    return single_result


class parallel_obj:
    """
    Empty class to hold parallel parameters.

    This is an internal class used by `parallel_routing` that never needs to be
    defined outside of that function.
    """

    def __init__(self):
        """Initialize attributes."""
        parallel_obj.particles = None
        parallel_obj.num_iter = None
        parallel_obj.Np_tracer = None
        parallel_obj.seed_xloc = None
        parallel_obj.seed_yloc = None


def parallel_routing(particles, num_iter, Np_tracer, seed_xloc, seed_yloc,
                     num_cores):
    """Do the parallel routing of particles.

    Function to do parallel routing of particles. Does this by duplicating the
    call to *Particle.run_iteration()* across different processes.

    **Inputs** :

        particle : :obj:`dorado.particle_track.Particles`
            An initialized :obj:`particle_track.Particles` object with some
            generated particles.

        num_iter : `int`
            Number of iterations to route particles for

        Np_tracer : `int`
            Number of particles to generate.

        seed_xloc : `list`
            List of x-coordinates over which to initially distribute the
            particles.

        seed_yloc : `list`
            List of y-coordinates over which to initially distribute the
            particles.

        num_cores : `int`
            Number of processors/cores to use for parallel part

    **Outputs** :

        par_result : `list`
            List of length(num_cores) with a dictionary of the beg/end indices
            and travel times for each particle computed by that process/core

    """
    # make parallel object to assign to function
    pobj = parallel_obj
    pobj.particles = particles
    pobj.num_iter = num_iter
    pobj.Np_tracer = Np_tracer
    pobj.seed_xloc = seed_xloc
    pobj.seed_yloc = seed_yloc

    # define list to pass to actual function that will be mapped to in parallel
    p_list = [pobj] * num_cores

    # create the parallel pool and run the process
    p = Pool(processes=num_cores)
    par_result = p.map(run_iter, p_list)
    p.terminate()

    return par_result
