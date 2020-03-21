# -*- coding: utf-8 -*-
"""
Particle tools to manage the internal functions related to the routing.

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

class Tools():



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

        # if sum of weights is above 0 normalize by sum of weights
        if np.nansum(weight_sfc) > 0:
            weight_sfc = weight_sfc / np.nansum(weight_sfc)
        # if sum of weight is above 0 normalize by sum of weights
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
        old_vel = self.velocity[ind[0],ind[1]]
        # new position velocity value
        new_vel = self.velocity[new_ind[0],new_ind[1]]
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
