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
        Particles class. e.g. particle = Particles(params)

        This initialization tries to assign each value from the parameter class,
        otherwise an error is raised or default values are assigned when
        possible/sensible

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
            raise ValueError("Stage array not specified")

        ### check if hydrodynamic model input has been specified
        if hasattr(params, 'model'):
            pass
        else:
            params.model = []

        ### Define x-component of discharge for all cells in domain
        try:
            if params.model == 'DeltaRCM':
                self.qx = params.qx
                self.qx[np.isnan(self.qx)] = 0
            else:
                self.qx = -1*params.qy
                self.qx[np.isnan(self.qx)] = 0
        except:
            raise ValueError("x-components of discharge values not specified")

        ### Define y-component of discharge for all cells in domain
        try:
            if params.model == 'DeltaRCM':
                self.qy = params.qy
                self.qy[np.isnan(self.qy)] = 0
            else:
                self.qy = params.qx
                self.qy[np.isnan(self.qy)] = 0
        except:
            raise ValueError("y-components of discharge values not specified")

        ### Define velocity field (for travel time calculation)
        # create modified depth array without 0s or nan values for division
        mod_depth = self.depth.copy()
        mod_depth[mod_depth==0] = 1e-6
        mod_depth[np.isnan(mod_depth)] = 1e-6
        # back out velocity field from discharge and depths
        # self.velocity = np.sqrt(self.qx**2+self.qy**2)*mod_depth/(mod_depth**2 + 1e-6)
        self.velocity = np.sqrt(self.qx**2+self.qy**2)/mod_depth
        # cannot have 0/nans - leads to infinite/nantravel times
        self.velocity[self.velocity==0] = 1e-6
        self.velocity[np.isnan(self.velocity)] = 1e-6


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
    def run_iteration(self,
                      start_xindices=None,
                      start_yindices=None,
                      start_times=None,
                      previous_walk_data=None,
                      target_time=None):
        '''
        Runs an iteration of the particle routing.
        Returns at each step the particle's locations and travel times.

        **Inputs** :

            start_xindices : `list`
                List of x locations to seed the particles [x1, x2, x3, ..., xn]
                if undefined, uses starting locations as given by the Particles
                class (seed_xloc)

            start_yindices : `list`
                List of y locations to seed the particles [y1, y2, y3, ..., yn]
                if undefined, uses starting locations as given by the Particles
                class (seed_yloc)

            start_times : `list`
                List of particle travel times [t1, t2, t3, ..., tn] if
                undefined, assumes no particles have travelled yet, so assigns
                zeros

            previous_walk_data : `list`
                Nested list of all prior x locations, y locations, and travel
                times (in that order). Order of indices is
                previous_walk_data[field][particle][iter], where e.g.
                [2][5][10] is the travel time of the 5th particle at the 10th
                iteration

            target_time : `float`
                The travel time (seconds) each particle should aim to have at
                end of this iteration. If left undefined, then just one
                iteration is run and the particles will be out of sync in time

        **Outputs** :

            all_walk_data : `list`
                Nested list of all x and y locations and travel times, with
                details same as input previous_walk_data

        '''

        if(previous_walk_data is not None):
            # If particle tracking has been run before, feed previous output array back into input
            # If this array exists, it overrides any starting indices given in function call
            all_xinds = previous_walk_data[0] # all previous locations
            start_xindices = [all_xinds[i][-1] for i in range(self.Np_tracer)] # most recent locations
            all_yinds = previous_walk_data[1]
            start_yindices = [all_yinds[i][-1] for i in range(self.Np_tracer)]
            all_times = previous_walk_data[2]
            start_times = [all_times[i][-1] for i in range(self.Np_tracer)]
        else:
            # if start locations not defined, then randomly assign them
            if start_xindices == None:
                start_xindices = map(lambda x: self.random_pick_seed(self.seed_xloc),
                                            range(self.Np_tracer)) # set starting x-index for all tracers
            if start_yindices == None:
                start_yindices = map(lambda x: self.random_pick_seed(self.seed_yloc),
                                            range(self.Np_tracer)) # set starting y-index for all tracers
            # initialize travel times list
            if start_times == None:
                start_times = [0.]*self.Np_tracer
            # Now initialize vectors that will create the structured list
            all_xinds = [[start_xindices[i]] for i in range(self.Np_tracer)]
            all_yinds = [[start_yindices[i]] for i in range(self.Np_tracer)]
            all_times = [[start_times[i]] for i in range(self.Np_tracer)]

		# If particles were placed near inlet and are having trouble starting motion, uncomment this:
        # self.qxn.flat[start_xindices] += 1 # add 1 to x-component of discharge at the start location
        # self.qyn.flat[start_yindices] += 1 # add 1 to y-component of discharge at the start location

        # merge x and y indices into list of [x,y] pairs
        start_pairs = [[start_xindices[i], start_yindices[i]] for i in range(self.Np_tracer)]

        # Do the particle movement
        if target_time == None:
            # If we're not aiming for a specific time, run a single iteration
            new_inds, travel_times = self.single_iteration(start_pairs, start_times)

            for ii in range(self.Np_tracer):
                all_xinds[ii].append(new_inds[ii][0]) # Append new information
                all_yinds[ii].append(new_inds[ii][1])
                all_times[ii].append(travel_times[ii])

            all_walk_data = [all_xinds, all_yinds, all_times] # Store travel information

        else: # If we ARE aiming for a specific time, iterate each particle until we get there
            # Loop through all particles
            for ii in range(self.Np_tracer):
                if(previous_walk_data is not None):
                    est_next_dt = all_times[ii][-1] - all_times[ii][-2]
                else:
                    est_next_dt = 0.1 # Initialize a guess for the next iteration's timestep
                count = 1

                # Loop until |target time - current time| < |target time - estimated next time|
                while abs(all_times[ii][-1] - target_time) >= abs(all_times[ii][-1] + est_next_dt - target_time):
                    # for particle ii, take a step from most recent index/time
                    new_inds, travel_times = self.single_iteration([[all_xinds[ii][-1], all_yinds[ii][-1]]],
                                                                   [all_times[ii][-1]])
                    all_xinds[ii].append(new_inds[0][0])
                    all_yinds[ii].append(new_inds[0][1])
                    all_times[ii].append(travel_times[0])

                    # Use that timestep to estimate how long the next one will take
                    est_next_dt = max(0.1, all_times[ii][-1] - all_times[ii][-2])
                    count += 1
                    if count > 1e4:
                        print('Warning: Particle iterations exceeded limit before reaching target time. Try smaller time-step')
                        break

            all_walk_data = [all_xinds, all_yinds, all_times]

        return all_walk_data
