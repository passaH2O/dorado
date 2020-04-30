# -*- coding: utf-8 -*-
"""
Parallel functionality using multiprocessing (included with Python).
For local parallelization using CPU cores.

Project Homepage: https://github.com/
"""

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

from particlerouting.particle_track import Particle
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

        data : `dict`
            Dictionary with beginning/end indices for particles and their
            associated travel times

    '''
    data = dict()
    particle = Particle(params)
    # do iterations
    for i in range(0,params.num_iter):
        if i == 0:
            start_inds, end_inds, travel_times = particle.run_iteration()
            beg_inds = start_inds # keep names consistent
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

    data['beg_inds'] = beg_inds
    data['end_inds'] = end_inds
    data['travel_times'] = travel_times

    return data



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

        single_result : `dict`
            Single dictionary with beg/end indices and travel times for all the
            particles

    '''

    # initiate final results dictionary
    single_result = dict()
    single_result['beg_inds'] = []
    single_result['end_inds'] = []
    single_result['travel_times'] = []

    # populate the dictionary
    for i in range(0,len(par_result)):
        if i == 0:
            single_result['beg_inds'] = par_result[i]['beg_inds']
            for j in range(0,len(par_result)):
                single_result['end_inds'].append(list(par_result[i]['end_inds'][j]))
            single_result['travel_times'] = par_result[i]['travel_times']
        else:
            for j in range(0,len(result[i]['beg_inds'])):
                single_result['beg_inds'].append([par_result[i]['beg_inds'][j][0],par_result[i]['beg_inds'][j][1]])
                single_result['end_inds'].append([par_result[i]['end_inds'][j][0],par_result[i]['end_inds'][j][1]])
            single_result['travel_times'] = single_result['travel_times'] + par_result[i]['travel_times']

    return single_result



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
