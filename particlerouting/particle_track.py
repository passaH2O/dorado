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
import time
from .particle_tools import Tools

class Particle(Tools):
    # class initialization, automatically run when class is defined
    # takes in attributes from a params object passed to it
    # e.g. testparticle = Particle(params)
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

        ### Define x-component of discharge for all cells in domain
        try:
            self.qx = params.qx
        except:
            raise ValueError("x-components of discharge values not specified")

        ### Define y-component of discharge for all cells in domain
        try:
            self.qy = params.qy
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

        ### Number of iterations for parcel routing when run_iteration is called
        try:
            self.itmax = params.itmax
        except:
            print('Max iterations not specified - using 1')
            self.itmax = 1 # max number of iterations is 1 if undefined

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
    def run_iteration(self,start_xindices=None,start_yindices=None,start_times=None):
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

        Outputs :
                    start_pairs : list [], of [x,y] pairs of the particle
                                  locations at the beginning of the iteration

                    new_inds : list [], of the new [x,y] locations for all of
                               the particles at the end of the iteration

                    travel_times : list [], of the travel times for each
                                   particle at the end of the timestep
        '''

        iter = 0 # set iteration counter to 0

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

        # loop until iterations are completed
        while (np.sum(current_inds) > 0) & (iter < self.itmax):

            iter += 1 # add +1 to the iter counter

            inds = current_inds #np.unravel_index(current_inds, self.depth.shape) # get indices as coordinates in the domain
            inds_tuple = [(inds[i][0], inds[i][1]) for i in range(len(inds))] # split the indices into tuples

            new_cells = map(lambda x: self.get_weight(x)
                            if x != (0,0) else 4, inds_tuple) # for each particle index get the weights

            new_inds = map(lambda x,y: self.calculate_new_ind(x,y)
                            if y != 4 else 0, inds_tuple, new_cells) # for each particle get the new index

            dist = map(lambda x,y,z: self.step_update(x,y,z) if x > 0
                       else 0, current_inds, new_inds, new_cells) # move each particle to the new index

            new_inds = np.array(new_inds, dtype = np.int) # put new indices into array
            new_inds[np.array(dist) == 0] = 0

            new_inds = self.check_for_boundary(new_inds,inds) # see if the indices are at boundaries

            # add the travel times
            temp_travel = map(lambda x,y: self.calc_travel_times(x,y) if x > 0
                                else 0, current_inds, new_inds)
            travel_times = [travel_times[i] + temp_travel[i] for i in range(0,len(travel_times))] # add to existing times
            travel_times = list(travel_times)

            # update current inds to the new ones
            current_inds = new_inds.tolist()

        return start_pairs, new_inds, travel_times
