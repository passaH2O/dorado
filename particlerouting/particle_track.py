# -*- coding: utf-8 -*-
"""
Particle class for managing the definition of particle attributes and parameters
of the domain as well as iterative movement of the particles through the domain.

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
import copy
import time
from .particle_tools import Tools


class params:
    '''
    The parameters class, `params`, which must be populated with user-defined
    attributes of the grid the particles will be modeled on.


    Expected Parameters:
    --------------------

        seed_xloc (`list`) : List of x-coordinates over which to initially
                             distribute the particles

        seed_yloc (`list`) : List of y-coordinates over which to initially
                             distribute the particles

        Np_tracer (`int`) : Number of particles to use

        dx (`float`) : Length along one square cell face

        depth (`numpy.ndarray`) : Array of water depth values, if absent then
                                  the stage and topography arrays will be used
                                  to compute it

        stage (`numpy.ndarray`) : Array of water stage values, if absent then
                                  the depth and topography arrays will be used
                                  to compute it

        qx (`numpy.ndarray`) : Array of the x-component of flow discharge

        qy (`numpy.ndarray`) : Array of the y-component of flow discharge

        u (`numpy.ndarray`) : Array of the x-component of flow velocity

        v (`numpy.ndarray`) : Array of the y-component of flow velocity


    Optional Parameters:
    --------------------

        topography (`numpy.ndarray`) : Array of cell elevation values

        model (`str`) : Name of the hydrodynamic model input being used
                        (e.g. 'DeltaRCM')

        theta (`float`) : First of two weighting parameters for the weighted
                          random walk. Default value is 1.0, higher values
                          give higher weighting probabilities to cells with
                          greater water depths

        gamma (`float`) : Second of two weighting parameters for the weighted
                          random walk. Default value is 0.05. Gamma must be in
                          the range [0,1]. Gamma == 0 means that the random
                          walk weights are independent of the discharge values,
                          and instead are based on the water surface gradient
                          (the stage). Gamma == 1 means that the random walk
                          weights are not dependent on the surface gradient,
                          and instead are based on the inertial forces (the
                          flow discharge).

        dry_depth (`float`) : Minimum depth for a cell to be considered wet,
                              default value is 0.1m

        cell_type (`numpy.ndarray`) : Array of the different types of cells in
                                      the domain where 2 = land, 1 = channel, 0
                                      = ocean, and -1 = edge. If not explicitly
                                      defined then the values are estimated
                                      based on the depth array and the defined
                                      dry_depth

        steepest_descent (`bool`) : Toggle for routing based on a steepest
                                    descent rather than the weighted random
                                    walk. If True, then the highest weighted
                                    cells are used to route the particles.
                                    Default value is False.


    This list of expected parameter values can also be obtained by querying the
    class attributes with `params.__dict__` or `vars(params)`

    '''

    def __init__(self):
        '''
        Creation of the expected variables for the params class.

        Variables are initialized as NoneType they need to be assigned by the
        user. Due to the wide variety of data formats produced by differeny
        hydrodynamic models, this is process is not automated and must be
        handled on a case-by-case basis.
        '''

        self.seed_xloc = None
        self.seed_yloc = None
        self.Np_tracer = None
        self.dx = None
        self.depth = None
        self.stage = None
        self.qx = None
        self.qy = None
        self.u = None
        self.v = None


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
        if getattr(params,'seed_xloc',None) is None:
            raise ValueError("No tracer seeding x-locations (params.seed_xloc)"
                             " have been defined")
        else:
            self.seed_xloc = params.seed_xloc

        if getattr(params,'seed_yloc',None) is None:
            raise ValueError("No tracer seeding y-locations (params.seed_yloc)"
                             " have been defined")
        else:
            self.seed_yloc = params.seed_yloc

        ### Define the number of tracers to be simulated
        if getattr(params,'Np_tracer',None) is None:
            raise ValueError("Number of tracer particles (params.Np_tracer)"
                             " has not been defined")
        else:
            self.Np_tracer = params.Np_tracer

        ### Define the length along one cell face (assuming square cells)
        if getattr(params,'dx',None) is None:
            raise ValueError("Length of cell face (params.dx) is undefined")
        else:
            self.dx = params.dx

        ### Define the water depth array
        if getattr(params,'depth',None) is not None:
            try:
                self.depth = params.depth
                self.depth[np.isnan(self.depth)] = 0
            except:
                raise ValueError("Water depth array incorrectly defined.")
        elif getattr(params,'stage',None) is not None and getattr(params,'topography',None) is not None:
            try:
                self.depth = params.stage - params.topography
                self.depth[self.depth < 0] = 0
                self.depth[np.isnan(self.depth)] = 0
            except:
                raise ValueError("Insufficient information: Specify depth")
        else:
            raise ValueError("Insufficient information: Specify depth")

        ### Define the water stage array
        if getattr(params,'stage',None) is not None:
            try:
                self.stage = params.stage
                self.stage[self.depth == 0] = np.nan
            except:
                raise ValueError("Water stage array incorrectly defined.")
        elif getattr(params,'topography',None) is not None and getattr(params,'depth',None) is not None:
            try:
                self.stage = params.topography + params.depth
                self.stage[self.depth == 0] = np.nan
            except:
                raise ValueError("Insufficient information: Specify stage")
        else:
            raise ValueError("Insufficient information: Specify stage")

        ### check if hydrodynamic model input has been specified
        if getattr(params, 'model', None) != None:
            pass
        else:
            params.model = []

        ### Define discharge and velocities for all cells in domain
        if params.model == 'DeltaRCM':
            if params.qx is not None and params.qy is not None:
                self.qx = params.qx
                self.qx[np.isnan(self.qx)] = 0
                self.u = self.qx*self.depth/(self.depth**2 + 1e-6)

                self.qy = params.qy
                self.qy[np.isnan(self.qy)] = 0
                self.v = self.qy*self.depth/(self.depth**2 + 1e-6)

            elif params.u is not None and params.v is not None:
                self.u = params.u
                self.u[np.isnan(self.u)] = 0
                self.qx = self.u*self.depth

                self.v = params.v
                self.v[np.isnan(self.v)] = 0
                self.qy = self.v*self.depth

            else:
                raise ValueError("Insufficient information: Specify velocities/discharge")
        else:
            if params.qx is not None and params.qy is not None:
                self.qx = -1*params.qy
                self.qx[np.isnan(self.qx)] = 0
                self.u = self.qx*self.depth/(self.depth**2 + 1e-6)

                self.qy = params.qx
                self.qy[np.isnan(self.qy)] = 0
                self.v = self.qy*self.depth/(self.depth**2 + 1e-6)

            elif params.u is not None and params.v is not None:
                self.u = -1*params.v
                self.u[np.isnan(self.u)] = 0
                self.qx = self.u*self.depth

                self.v = params.u
                self.v[np.isnan(self.v)] = 0
                self.qy = self.v*self.depth

            else:
                raise ValueError("Insufficient information: Specify velocities/discharge")

        ### Define field of velocity magnitude (for travel time calculation)
        self.velocity = np.sqrt(self.u**2+self.v**2)
        # cannot have 0/nans - leads to infinite/nantravel times
        self.velocity[self.velocity < 1e-6] = 1e-6


        ########## OPTIONAL PARAMETERS (Have default values) ##########
        ### Define the theta used to weight the random walk
        # Higher values give higher weighting probabilities to deeper cells
        try:
            self.theta = params.theta
        except:
            print("Theta parameter not specified - using 1.0")
            self.theta = 1.0 # if unspecified use 1

        ### Gamma parameter used to weight the random walk
        # Sets weight ratio (between 0 and 1):
                            # 1 = water surface gradient only (stage based)
                            # 0 = inertial force only (discharge based)
        try:
            self.gamma = params.gamma
        except:
            print("Gamma parameter not specified - using 0.05")
            self.gamma = 0.05

        ### Minimum depth for cell to be considered wet
        try:
            self.dry_depth = params.dry_depth
        except:
            print("minimum depth for wetness not defined - using 10 cm")
            self.dry_depth = 0.1

        ### Cell types: 2 = land, 1 = channel, 0 = ocean, -1 = edge
        try:
            self.cell_type = params.cell_type
        except:
            print("Cell Types not specified - Estimating from depth")
            self.cell_type = np.zeros_like(self.depth, dtype='int')
            self.cell_type[self.depth < self.dry_depth] = 2
            self.cell_type = np.pad(self.cell_type[1:-1,1:-1], 1, 'constant', constant_values = -1)

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
                iteration is run and the particles will be out of sync in time.
                Note that this loop will terminate before the target_time if
                the particle exceeds the hard-coded limit of 1e4 steps

        **Outputs** :

            all_walk_data : `list`
                Nested list of all x and y locations and travel times, with
                details same as input previous_walk_data

        '''

        if previous_walk_data is not None:
            # If particle tracking has been run before, feed previous output array back into input
            # If this array exists, it overrides any starting indices given in function call
            all_xinds = previous_walk_data[0] # all previous locations
            start_xindices = [all_xinds[i][-1] for i in list(range(self.Np_tracer))] # most recent locations
            all_yinds = previous_walk_data[1]
            start_yindices = [all_yinds[i][-1] for i in list(range(self.Np_tracer))]
            all_times = previous_walk_data[2]
            start_times = [all_times[i][-1] for i in list(range(self.Np_tracer))]
        else:
            # if start locations not defined, then randomly assign them
            if start_xindices == None:
                start_xindices = [self.random_pick_seed(self.seed_xloc) for x in list(range(self.Np_tracer))] # set starting x-index for all tracers
            if start_yindices == None:
                start_yindices = [self.random_pick_seed(self.seed_yloc) for x in list(range(self.Np_tracer))] # set starting y-index for all tracers
            # initialize travel times list
            if start_times == None:
                start_times = [0.]*self.Np_tracer
            # Now initialize vectors that will create the structured list
            all_xinds = [[start_xindices[i]] for i in list(range(self.Np_tracer))]
            all_yinds = [[start_yindices[i]] for i in list(range(self.Np_tracer))]
            all_times = [[start_times[i]] for i in list(range(self.Np_tracer))]

        # merge x and y indices into list of [x,y] pairs
        start_pairs = [[start_xindices[i], start_yindices[i]] for i in list(range(self.Np_tracer))]

        # Do the particle movement
        if target_time == None:
            # If we're not aiming for a specific time, run a single iteration
            new_inds, travel_times = self.single_iteration(start_pairs, start_times)

            for ii in list(range(self.Np_tracer)):
                # Don't duplicate location if particle is standing still at a boundary
                if new_inds[ii] != start_pairs[ii]:
                    all_xinds[ii].append(new_inds[ii][0]) # Append new information
                    all_yinds[ii].append(new_inds[ii][1])
                    all_times[ii].append(travel_times[ii])

            all_walk_data = [all_xinds, all_yinds, all_times] # Store travel information

        else: # If we ARE aiming for a specific time, iterate each particle until we get there
            # Loop through all particles
            for ii in list(range(self.Np_tracer)):
                if previous_walk_data is not None:
                    est_next_dt = all_times[ii][-1] - all_times[ii][-2]
                else:
                    est_next_dt = 0.1 # Initialize a guess for the next iteration's timestep
                count = 1

                # Only iterate if this particle isn't already at a boundary:
                if -1 not in self.cell_type[all_xinds[ii][-1]-1:all_xinds[ii][-1]+2,
                                            all_yinds[ii][-1]-1:all_yinds[ii][-1]+2]:
                    # Loop until |target time - current time| < |target time - estimated next time|
                    while abs(all_times[ii][-1] - target_time) >= abs(all_times[ii][-1] + est_next_dt - target_time):
                        # for particle ii, take a step from most recent index/time
                        new_inds, travel_times = self.single_iteration([[all_xinds[ii][-1], all_yinds[ii][-1]]],
                                                                       [all_times[ii][-1]])

                        # Don't duplicate location if particle is standing still at a boundary
                        if new_inds[0] != [all_xinds[ii][-1], all_yinds[ii][-1]]:
                            all_xinds[ii].append(new_inds[0][0])
                            all_yinds[ii].append(new_inds[0][1])
                            all_times[ii].append(travel_times[0])
                        else:
                            break

                        # Use that timestep to estimate how long the next one will take
                        est_next_dt = max(0.1, all_times[ii][-1] - all_times[ii][-2])
                        count += 1
                        if count > 1e4:
                            print('Warning: Particle iterations exceeded limit before reaching target time. Try smaller time-step')
                            break

            all_walk_data = [all_xinds, all_yinds, all_times]

        return all_walk_data
