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

    def __init__(self, params):
        '''
        Methods require a set of parameters to be assigned
        Try to assign each value from the parameter file, otherwise raise error
        '''
        try:
            self.inlet = params.inlet
        except:
            raise ValueError("No inlet location defined")

        try:
            self.Np_water = params.Np_water
        except:
            raise ValueError("Number of water parcels not specified")

        try:
            self.Qp_water = params.Qp_water
        except:
            raise ValueError("Volume of each water parcel not specified")

        try:
            self.dx = params.dx
        except:
            raise ValueError("Cell size not specified")

        try:
            self.qxn = params.qxn
        except: # might be able to define w/o being specified
            raise ValueError("x components of discharge not specified")

        try:
            self.qyn = params.qyn
        except: # might be able to define w/o being specified
            raise ValueError("y components of discharge not specified")

        try:
            self.qwn = params.qwn
        except: # might be able to define w/o being specified
            raise ValueError("Per parcel discharge values not specified???")

        try:
            self.indices = params.indices
        except:
            # should be able to write method to define empty array for indices
            pass

        try:
            self.looped = params.looped
        except:
            # should be able to write method to define empty array if not given
            pass

        try:
            self.itmax = params.itmax
        except:
            self.itmax = 1 # max number of iterations is 1 if undefined

        try:
            self.depth = params.depth
        except:
            raise ValueError("Depth array not specified")

        try:
            self.free_surf_flag = params.free_surf_flag
        except:
            raise ValueError("still dont know what this is ..")

        try:
            self.stage = params.stage
        except:
            raise ValueError("Stage values not specified")

        try:
            self.distances = params.distances
        except:
            # can define this if not given
            pass

        try:
            self.qx = params.qx
        except:
            raise ValueError("x-components of discharge values not specified")

        try:
            self.qy = params.qy
        except:
            raise ValueError("y-components of discharge values not specified")

        try:
            self.ivec = params.ivec
        except:
            # can define this if not given
            pass

        try:
            self.jvec = params.jvec
        except:
            # can define this if not given
            pass

        try:
            self.dry_depth = params.dry_depth
        except:
            raise ValueError("minimum depth for wetness not defined")

        try:
            self.gamma = params.gamma
        except:
            raise ValueError("parameter gamma not specified")

        try:
            self.iwalk = params.iwalk
        except:
            # can be defined if not given
            pass

        try:
            self.jwalk = params.jwalk
        except:
            # can be defined if not given
            pass

        try:
            self.L = params.L
        except:
            self.L = np.shape(qx,1) # check syntax

        try:
            self.W = params.W
        except:
            self.W = np.shape(qx,2) # check syntax

        try:
            self.cell_type = params.cell_type
        except:
            raise ValueError("Cell Types not specified")

        try:
            self.L0 = params.L0
        except:
            # might be able to get around this
            pass

        try:
            self.CTR = params.CTR
        except:
            self.CTR = np.round(self.W/2)

        try:
            self.sfc_visit = params.sfc_visit
        except:
            self.sfc_visit = np.zeros_like(self.depth)

        try:
            self.sfc_sum = params.sfc_sum
        except:
            self.sfc_sum = np.zeros_like(self.depth)

        try:
            self.theta_water = params.theta_water
        except:
            raise ValueError("Theta Water not defined")



    def init_water_iteration(self):

        self.qxn[:] = 0; self.qyn[:] = 0; self.qwn[:] = 0

        self.free_surf_flag[:] = 0
        self.indices[:] = 0
        self.sfc_visit[:] = 0
        self.sfc_sum[:] = 0

        self.pad_stage = np.pad(self.stage, 1, 'edge')

        self.pad_depth = np.pad(self.depth, 1, 'edge')

        self.pad_cell_type = np.pad(self.cell_type, 1, 'edge')



    def run_water_iteration(self):

        iter = 0 # set iteration counter to 0
        start_indices = map(lambda x: self.random_pick_inlet(self.inlet),
                                      range(self.Np_water)) # set starting index for the water parcel

        self.qxn.flat[start_indices] += 1 # add 1 to x-component of discharge at the start location
        self.qwn.flat[start_indices] += self.Qp_water / self.dx / 2. # initial per parcel water discharge

        self.indices[:,0] = start_indices # save the start index for the parcel
        current_inds = list(start_indices) # list of the start location indices

        self.looped[:] = 0 # set this loop toggle to 0


        while (sum(current_inds) > 0) & (iter < self.itmax):

            iter += 1 # add +1 to the iter counter

            inds = np.unravel_index(current_inds, self.depth.shape) # get indices as coordinates in the domain
            inds_tuple = [(inds[0][i], inds[1][i]) for i in range(len(inds[0]))] # split the indices into tuples


            new_cells = map(lambda x: self.get_weight(x)
                            if x != (0,0) else 4, inds_tuple) # for each particle index get the weights


            new_inds = map(lambda x,y: self.calculate_new_ind(x,y)
                            if y != 4 else 0, inds_tuple, new_cells) # for each particle get the new index


            dist = map(lambda x,y,z: self.step_update(x,y,z) if x > 0
                       else 0, current_inds, new_inds, new_cells) # move each particle to the new index


            new_inds = np.array(new_inds, dtype = np.int) # put new indices into array
            new_inds[np.array(dist) == 0] = 0


            self.indices[:,iter] = new_inds # assign the new index to indices array

            current_inds = self.check_for_loops(new_inds, iter) # run check for loops function

            current_inds = self.check_for_boundary(current_inds) # see if the indices are at boundaries

            self.indices[:,iter] = current_inds # assign indices as the current_inds list

            current_inds[self.free_surf_flag > 0] = 0 # check this free surface flag ???




    ### random pick inlet cell function (can be adapted to control particle drop location)
    def random_pick_inlet(self, choices, probs = None):
        '''
        Randomly pick a number from array choices weighted by array probs
        Values in choices are column indices
        Return a tuple of the randomly picked index for row 0
        '''
        # randomly pick the inlet cell to use given a list of inlet indices
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
        # define weights as least 0 or a greater value
        weight_sfc = np.maximum(0,
                     (self.stage[ind] - stage_ind) / self.distances)
        # define weight integers as at least 0 or greater value
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

        # define actual weight by using gamma, and the defined weights
        self.weight = self.gamma * weight_sfc + (1 - self.gamma) * weight_int
        # modify the weight by the depth and theta_water parameter
        self.weight = depth_ind ** self.theta_water * self.weight
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
        new_ind_flat = np.ravel_multi_index(new_ind, self.depth.shape)

        return new_ind_flat




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
            self.qxn.flat[ind] += jstep / dist
            self.qyn.flat[ind] += istep / dist
            self.qwn.flat[ind] += self.Qp_water / self.dx / 2.
            # identify the new index for discharge
            # add the distance to the new values
            self.qxn.flat[new_ind] += jstep / dist
            self.qyn.flat[new_ind] += istep / dist
            self.qwn.flat[new_ind] += self.Qp_water / self.dx / 2.

        return dist




    ### check loops
    def check_for_loops(self, inds, it):

        looped = [len(i[i>0]) != len(set(i[i>0])) for i in self.indices]

        for n in range(self.Np_water):

            ind = inds[n]

            if (looped[n]) and (ind > 0) and (it > self.L0):

                self.looped[n] += 1

                it = np.unravel_index(ind, self.depth.shape)

                px, py = it

                Fx = px - 1
                Fy = py - self.CTR

                Fw = np.sqrt(Fx**2 + Fy**2)

                if Fw != 0:
                    px = px + np.round(Fx / Fw * 5.)
                    py = py + np.round(Fy / Fw * 5.)

                px = max(px, self.L0)
                px = int(min(self.L - 2, px))

                py = max(1, py)
                py = int(min(self.W - 2, py))

                nind = np.ravel_multi_index((px,py), self.depth.shape)

                inds[n] = nind

                self.free_surf_flag[n] = -1


        return inds




    ### bndy check
    def check_for_boundary(self, inds):
        # check if the cell is on an edge (type==-1) if it is on and edge AND it is
        # also conditional on a free_surf_flag??? not sure that that flag is
        self.free_surf_flag[(self.cell_type.flat[inds] == -1) & (self.free_surf_flag == 0)] = 1

        self.free_surf_flag[(self.cell_type.flat[inds] == -1) & (self.free_surf_flag == -1)] = 2

        inds[self.free_surf_flag == 2] = 0 # if free surface flag is 2, make that ind value 0

        return inds




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
