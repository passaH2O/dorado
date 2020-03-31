# -*- coding: utf-8 -*-
"""
Particle class for managing the definition of particle attributes and parameters
of the domain as well as iterative movement of the particles through the domain.

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
import copy
import time
from .particle_tools import Tools

class Particle(Tools):
    '''
    Class for the particle or set of particles that is going to be routed
    '''
    def __init__(self, params):
        '''
        Methods require a class of parameters (params) to be passed to the
        Particles class i.e. particle = Particles(params)

        Try to assign each value from the parameter file, otherwise raise error
        or default values are assigned when possible/sensible
        '''

        ########## REQUIRED PARAMETERS ##########
        ### Define the seeding locations as list of x and y locations
        try:
            self.seed_xloc = params.seed_xloc
        except:
            raise ValueError("No tracer seeding x-locations defined")

        try:
            self.seed_yloc = params.seed_yloc
        except:
            raise ValueError("No tracer seeding y-locations defined")

        ### Define the number of tracers to be simulated
        try:
            self.Np_tracer = params.Np_tracer
        except:
            raise ValueError("Number of tracer particles not specified")

        ### Define the length along one cell face (assuming square cells)
        try:
            self.dx = params.dx
        except:
            raise ValueError("Cell size not specified")

        ### Define the water depth array
        try:
            self.depth = params.depth
        except:
            raise ValueError("Depth array not specified")

        ### Define the water stage array
        try:
            self.stage = params.stage
        except:
            print("Stage values not specified - using depth values")
            self.stage = self.depth

        ### check if hydrodynamic model input has been specified
        if hasattr(params, 'model'):
            pass
        else:
            params.model = []

        ### Define x-component of discharge for all cells in domain
        try:
            if params.model == 'DeltaRCM':
                self.qx = params.qx
            else:
                self.qx = -1*params.qy
        except:
            raise ValueError("x-components of discharge values not specified")

        ### Define y-component of discharge for all cells in domain
        try:
            if params.model == 'DeltaRCM':
                self.qy = params.qy
            else:
                self.qy = params.qx
        except:
            raise ValueError("y-components of discharge values not specified")

        ### Define velocity field (for travel time calculation)
        # create modified depth array without 0s or nan values for division
        mod_depth = self.depth.copy()
        mod_depth[mod_depth==0] = 1e-10
        mod_depth[np.isnan(mod_depth)] = 1e-10
        # back out velocity field from discharge and depths
        self.velocity = np.sqrt(self.qx**2+self.qy**2)/mod_depth
        # cannot have 0/nans - leads to infinite/nantravel times
        self.velocity[self.velocity==0] = 1e-10
        self.velocity[np.isnan(self.velocity)] = 1e-10


        ########## OPTIONAL PARAMETERS (Have default values) ##########
        ### Define the theta used to weight the random walk
        try:
            self.theta = params.theta
        except:
            print('Theta not specified - using 1.0')
            self.theta = 1.0 # if unspecified use 1

        ### Minimum depth for cell to be considered wet
        try:
            self.dry_depth = params.dry_depth
        except:
            print("minimum depth for wetness not defined - using 0.1")
            self.dry_depth = 0.1

        ### Gamma parameter
        # sets weight ratio:
                            # either water surface slope (depth based)
                            # or
                            # inertial force (discharge based)
        try:
            self.gamma = params.gamma
        except:
            print("parameter gamma not specified - using 0.05")
            self.gamma = 0.05

        ### Cell types: 2 = land, 1 = channel, 0 = ocean, -1 = edge
        try:
            self.cell_type = params.cell_type
        except:
            print("Cell Types not specified - using zeros")
            self.cell_type = np.zeros(np.shape(self.stage))

        ### Steepest descent toggle - turns off randomness and uses highest
        # weighted value instead of doing weighted random walk
        # note: chooses randomly in event of ties
        try:
            if params.steepest_descent == True:
                print("Using steepest descent")
                self.steepest_descent = True
            else:
                print("Using weighted random walk")
                self.steepest_descent = False
        except:
            print("Using weighted random walk")
            self.steepest_descent = False


        ########## DEFAULT PARAMETERS (Can be defined otherwise) ##########

        sqrt2 = np.sqrt(2)
        sqrt05 = np.sqrt(0.5)

        ### Define distances between cells in D8 sense
        try:
            self.distances = params.distances
        except:
            # defined this if not given
            self.distances = np.array([[sqrt2, 1, sqrt2],
                                       [1, 1, 1],
                                       [sqrt2, 1, sqrt2]])

        ### D8 components of x-unit vector
        try:
            self.ivec = params.ivec
        except:
            # define this if not given
            self.ivec = np.array([[-sqrt05, 0, sqrt05],
                                  [-1, 0, 1],
                                  [-sqrt05, 0, sqrt05]])

        ### D8 components of y-unit vector
        try:
            self.jvec = params.jvec
        except:
            # define this if not given
            self.jvec = np.array([[-sqrt05, -1, -sqrt05],
                                  [0, 0, 0],
                                  [sqrt05, 1, sqrt05]])

        ### Positive/Negative x-directions
        try:
            self.iwalk = params.iwalk
        except:
            # defined if not given
            self.iwalk = np.array([[-1, 0, 1],
                                   [-1, 0, 1],
                                   [-1, 0, 1]])

        ### Positive/Negative y-directions
        try:
            self.jwalk = params.jwalk
        except:
            # defined if not given
            self.jwalk = np.array([[-1, -1, -1],
                                   [0, 0, 0],
                                   [1, 1, 1]])

        # establish some zero matrices that are filled in when iteration is run
        self.qxn = np.zeros(np.shape(self.stage))
        self.qyn = np.zeros(np.shape(self.stage))
        self.sfc_visit = np.zeros(np.shape(self.stage))
        self.sfc_sum = np.zeros(np.shape(self.stage))

        # pad stage and depth arrays to identify edges
        self.pad_stage = np.pad(self.stage, 1, 'edge')
        self.pad_depth = np.pad(self.depth, 1, 'edge')
        self.pad_cell_type = np.pad(self.cell_type, 1, 'constant', constant_values = -1)


    # run an iteration where particles are moved
    # have option of specifying the particle start locations
    # otherwise they are randomly placed within x and y seed locations
    def run_iteration(self,start_xindices=None,start_yindices=None,start_times=None,time_step=None):
        '''
        Runs an iteration of the particle routing.
        Returns the original particle locations and their final locations.

        Inputs :
                    start_xindices : list of x locations to seed the particles
                                     [x1, x2, x3, ..., xn]
                                     if undefined, uses starting locations as
                                     given by the Particles class (seed_xloc)

                    start_yindices : list of y locations to seed the particles
                                     [y1, y2, y3, ..., yn]
                                     if undefined, uses starting locations as
                                     given by the Particles class (seed_yloc)

                    start_times : list of particle travel times
                                  [t1, t2, t3, ..., tn]
                                  if undefined, assumes no particles have
                                  travelled yet, so assigns zeros

                    time_step : model timestep (seconds), or the travel time
                                each particle should have at end of iteration,
                                if left undefined, then just one iteration is
                                run and the particles will be out of sync
                                (time-wise)

        Outputs :
                    start_pairs : list [], of [x,y] pairs of the particle
                                  locations at the beginning of the iteration

                    new_inds : list [], of the new [x,y] locations for all of
                               the particles at the end of the iteration

                    travel_times : list [], of the travel times for each
                                   particle at the end of the timestep
        '''

        # if start locations not defined, then randomly assign them
        if start_xindices == None:
            start_xindices = map(lambda x: self.random_pick_seed(self.seed_xloc),
                                        range(self.Np_tracer)) # set starting x-index for all tracers
        if start_yindices == None:
            start_yindices = map(lambda x: self.random_pick_seed(self.seed_yloc),
                                        range(self.Np_tracer)) # set starting y-index for all tracers

        self.qxn.flat[start_xindices] += 1 # add 1 to x-component of discharge at the start location
        self.qyn.flat[start_yindices] += 1 # add 1 to y-component of discharge at the start location

        # merge x and y indices into list of [x,y] pairs
        start_pairs = [[start_xindices[i],start_yindices[i]] for i in range(0,len(start_xindices))]

        # copy list of start location indices so start_pairs can be kept unchanged
        current_inds = start_pairs

        # initialize travel times list
        if start_times == None:
            travel_times = np.zeros(len(current_inds))
        else:
            travel_times = start_times

        # do the particle movement
        if time_step == None:
            # run a single iteration
            new_inds, travel_times = self.single_iteration(current_inds, travel_times)

        else:
            # define end time for particles as average of their current time
            # plus the defined timestep duration
            end_time = np.mean(travel_times) + time_step

            # run a single particle iteration
            temp_inds, temp_times = self.single_iteration(current_inds, travel_times)

            # copy results from first iteration into final lists
            new_inds = copy.copy(temp_inds)
            travel_times = copy.copy(temp_times)

            # iterate more as long as smalles of the travel times is below
            # 90% of the end time
            while np.min(temp_times) < 0.9*end_time:
                temp_inds, temp_times = self.single_iteration(temp_inds, temp_times)
                # check each particle time
                for i in range(0,len(temp_times)):
                    # if the particle time is close to end time then save it
                    if temp_times[i] > 0.9*end_time:
                        new_inds[i] = temp_inds[i]
                        travel_times[i] = temp_times[i]


        return start_pairs, new_inds, travel_times
