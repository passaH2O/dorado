#! /usr/bin/env python
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


class Particle():
    # class initialization, automatically run when class is defined
    # takes in attributes from a params object passed to it\
    # e.g. testparticle = Particle(params)
    def __init__(self, params):
        '''
        Methods require a set of parameters to be assigned
        Try to assign each value from the parameter file, otherwise raise error
        or default values are assigned when possible
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
        self.velocity = np.sqrt(self.qx**2+self.qy**2)/self.depth
        # cannot have 0/nans - leads to infinite/nantravel times
        self.velocity[self.velocity==0] = 1e-10
        self.velocity[np.isnan(self.velocity)] = 1e-10

        ### Define the theta used to weight the random walk
        try:
            self.theta = params.theta
        except:
            raise ValueError("Theta weight not defined")


        ########## OPTIONAL PARAMETERS (Have default values) ##########
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
        Runs an iteration of the particle routing
        Returns the original particle locations and their final locations
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



    ### random pick seeding location
    def random_pick_seed(self, choices, probs = None):
        '''
        Randomly pick a number from array choices weighted by array probs
        Values in choices are column indices
        Return a tuple of the randomly picked index for row 0
        '''
        # randomly pick tracer drop cell to use given a list of potential spots
        if not probs:
            probs = np.array([1 for i in range(len(choices))])
        # find the corresponding index value from the input 'choices' list of indices
        cutoffs = np.cumsum(probs)
        idx = cutoffs.searchsorted(np.random.uniform(0, cutoffs[-1]))

        return choices[idx]



    ### pull weights
    def get_weight(self, ind):
        # pad stage array with 1's around it
        stage_ind = self.pad_stage[ind[0]-1+1:ind[0]+2+1, ind[1]-1+1:ind[1]+2+1]
        # define water surface gradient weight component (minimum of 0)
        weight_sfc = np.maximum(0,
                     (self.stage[ind] - stage_ind) / self.distances)
        # define flow inertial weighting component (minimum of 0)
        weight_int = np.maximum(0, (self.qx[ind] * self.jvec +
                                    self.qy[ind] * self.ivec) / self.distances)

        # if the value of the first index coord is 0, make weights 0
        # might prevent travel up/out of domain???
        if ind[0] == 0:
            weight_sfc[0,:] = 0
            weight_int[0,:] = 0

        # add 1 value of padding around
        depth_ind = self.pad_depth[ind[0]-1+1:ind[0]+2+1, ind[1]-1+1:ind[1]+2+1]
        ct_ind = self.pad_cell_type[ind[0]-1+1:ind[0]+2+1, ind[1]-1+1:ind[1]+2+1]
        # if the depth is below minimum depth for cell to be weight or it is a cell
        # type of 'land' which is the walls, then make it impossible for the parcel
        # to travel there by setting associated weight to 0
        weight_sfc[(depth_ind <= self.dry_depth) | (ct_ind == -2)] = 0
        weight_int[(depth_ind <= self.dry_depth) | (ct_ind == -2)] = 0

        # if sum of weights is above 0 normalize by sum?
        if np.nansum(weight_sfc) > 0:
            weight_sfc = weight_sfc / np.nansum(weight_sfc)
        # if sum of weight is above 0 normalize by sum?
        if np.nansum(weight_int) > 0:
            weight_int = weight_int / np.nansum(weight_int)

        # define actual weight by using gamma, and the defined weight components
        self.weight = self.gamma * weight_sfc + (1 - self.gamma) * weight_int
        # modify the weight by the depth and theta weighting parameter
        self.weight = depth_ind ** self.theta * self.weight
        # if the depth is below the minimum depth then location is not considered
        # therefore set the associated weight to nan
        self.weight[depth_ind <= self.dry_depth] = np.nan
        # randomly pick the new cell for the particle to move to using the
        # random_pick function and the set of weights just defined
        new_cell = self.random_pick(self.weight)

        return new_cell



    ### calculate new index
    def calculate_new_ind(self, ind, new_cell):
        # add the index and the flattened x and y walk component
        # x,y walk component is related to the next cell chosen as a 1-8 location
        new_ind = (ind[0] + self.jwalk.flat[new_cell], ind[1] +
                   self.iwalk.flat[new_cell])
        # using the new_ind value re-ravel the index into a properly shaped array
        #new_ind_flat = np.ravel_multi_index(new_ind, self.depth.shape)

        return new_ind#new_ind_flat



    ### update step
    def step_update(self, ind, new_ind, new_cell):
        # assign x-step by pulling 1-8 value from x-component walk 1-8 directions
        istep = self.iwalk.flat[new_cell]
        # assign y-step by pulling 1-8 value from y-component walk 1-8 directions
        jstep = self.jwalk.flat[new_cell]
        # compute the step distance to be taken
        dist = np.sqrt(istep**2 + jstep**2)
        # check to ensure the step is actually traversing some distance
        if dist > 0:
            # identify original discharge locations in flattened arrays
            # add the new location to the original index information
            self.qxn.flat[ind[0]] += jstep / dist
            self.qyn.flat[ind[1]] += istep / dist
            # identify the new index for discharge
            # add the distance to the new values
            self.qxn.flat[new_ind[0]] += jstep / dist
            self.qyn.flat[new_ind[1]] += istep / dist

        return dist



    ### calculate travel time using avg of velocity and old and new index
    def calc_travel_times(self, ind, new_ind):
        # get old position velocity value
        old_vel = self.velocity.flat[ind]
        # new position velocity value
        new_vel = self.velocity.flat[new_ind]
        # avg velocity
        avg_vel = np.mean([old_vel,new_vel])
        # travel time based on cell size and mean velocity
        trav_time = 1/avg_vel * self.dx

        return trav_time


    ### bndy check
    def check_for_boundary(self, new_inds, current_inds):
        # check if the new indices are on edges (type==-1)
        # if so then don't let parcel go there
        for i in range(0,len(new_inds)):
            if self.pad_cell_type[new_inds[i][0],new_inds[i][1]] == -1:
                new_inds[i][0] = current_inds[i][0]
                new_inds[i][1] = current_inds[i][1]
            else:
                pass

        return new_inds



    ### random pick from weighted array probabilities
    def random_pick(self, probs):
        '''
        Randomly pick a number weighted by array probs (len 8)
        Return the index of the selected weight in array probs
        '''
        # check for the number of nans in the length 8 array of locations around the location
        num_nans = sum(np.isnan(probs))
        # if there are no nans, then everywhere there is no nan in the probs list is assigned a 1
        if np.nansum(probs) == 0:
            probs[~np.isnan(probs)] = 1
            probs[1,1] = 0 # except location 1,1 which is assigned a 0

        probs[np.isnan(probs)] = 0 # any nans are assigned as 0
        cutoffs = np.cumsum(probs) # cumulative sum of all probabilities
        # randomly pick indices from cutoffs based on uniform distribution
        idx = cutoffs.searchsorted(np.random.uniform(0, cutoffs[-1]))

        return idx
