# -*- coding: utf-8 -*-
"""
Parallel functionality using multiprocessing (included with Python).
For local parallelization using CPU cores.

Project Homepage: https://github.com/
"""
from __future__ import division, print_function, absolute_import
from builtins import range, map
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

from .particle_track import Particle
from multiprocessing import Pool


# convert run script into a function that returns a single dictionary 'data'
def run_iter(params):
    '''
    Uses params class to define a *Particle* and does iterations of the
    particle. Requires a new params variable : *params.num_iter* to define the
    number of iterations to route the particles for.

    **Inputs** :

        params : `obj`
            Class with all the parameters for the particle routing

    **Outputs** :

        all_walk_data : `list`
            Nested list of all x and y locations and travel times, with
            details same as input previous_walk_data

    '''
    particle = Particle(params)
    all_walk_data = None
    # do iterations
    for i in list(range(0,params.num_iter)):
        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data)

    return all_walk_data



def combine_result(par_result):
    '''
    Take the parallel resulting list and combine the entries so that a single
    dictionary with the beginning/end indices and travel times for all particles
    is returned as opposed to a list with one dictionary per parallel process

    **Inputs** :

        par_result : `list`
            List of length(num_cores) with a dictionary of the beg/end indices
            and travel times for each particle computed by that process/core

    **Outputs** :

        single_result : `list`
            Nested list that matches 'all_walk_data'

    '''

    ## list version
    # initiate final results list [[xinds],[yinds],[times]]
    single_result = [[],[],[]]

    # populate the dictionary
    # loop through results for each core
    for i in list(range(0,len(par_result))):
        # append results for each category
        for j in list(range(0,len(par_result[i][0]))):
            single_result[0].append(par_result[i][0][j])
            single_result[1].append(par_result[i][1][j])
            single_result[2].append(par_result[i][2][j])


    return single_result

    ### dictionary version
    # # initiate final results dictionary
    # single_result = dict()
    # single_result['x_inds'] = []
    # single_result['y_inds'] = []
    # single_result['travel_times'] = []
    #
    # # populate the dictionary
    # # loop through results for each core
    # for i in range(0,len(par_result)):
    #     # append results for each category
    #     for j in range(0,len(par_result[i][0])):
    #         single_result['x_inds'].append(par_result[i][0][j])
    #         single_result['y_inds'].append(par_result[i][1][j])
    #         single_result['travel_times'].append(par_result[i][2][j])
    #
    #
    # return single_result



def parallel_routing(params,num_iter,num_cores):
    '''
    Function to do parallel routing of particles. Does this by duplicating the
    call to *Particle.run_iteration()* across different processes.

    **Inputs** :

        params : `obj`
            Normal parameters for the particles in a class

        num_iter : `int`
            Number of iterations to route particles for

        num_cores : `int`
            Number of processors/cores to use for parallel part

    **Outputs** :

        par_result : `list`
            List of length(num_cores) with a dictionary of the beg/end indices
            and travel times for each particle computed by that process/core

    '''

    # assign number of iterations to the params class
    params.num_iter = num_iter

    # make params list to apply to map function
    params_list = [params]
    params_list = params_list * num_cores

    # create the parallel pool and run the process
    p = Pool(processes=num_cores)
    par_result = p.map(run_iter, params_list)
    p.terminate()

    return par_result
