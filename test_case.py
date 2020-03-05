# make a test case to try the particle class

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

# define an empty class
class pobj():
    pass

# define square root variables
sqrt2 = np.sqrt(2)
sqrt05 = np.sqrt(0.5)

# create params and then assign the parameters
params = pobj()

# pull some from a deltarcm output so stage is varied
rcm_output = Dataset('jgrdup_50_pydeltarcm.nc','r')

# pull stage
rcm_stage = rcm_output.variables['stage'][:]
params.stage = rcm_stage.data[-1,:,:]
# plt.figure()
# plt.imshow(params.stage)
# plt.show()

# pull depth
rcm_depth = rcm_output.variables['depth'][:]
params.depth = rcm_depth.data[-1,:,:]
# plt.figure()
# plt.imshow(params.depth)
# plt.show()

params.seed_xloc = list(range(17,23))
params.seed_yloc = list(range(126,131))
params.Np_tracer = 250
params.dx = 50.
params.qx = np.zeros(np.shape(params.depth))
params.qy = np.zeros(np.shape(params.depth))
params.theta = 1.0
# params.gamma = 10.0
params.itmax = 50

# create some discharge -- doesn't matter? stage used to weight?
# params.qx[40:60,40:60] = 100.0
# params.qy[30:45,30:45] = 0.0
# plt.figure()
# plt.imshow(np.sqrt(params.qx**2+params.qy**2))
# plt.colorbar()
# plt.title('Discharge field')
# plt.show()

# try running it
from particle_track import Particle

test = Particle(params)

test.init_water_iteration()

# # make pre-iteration plot
# plt.figure()
# plt.subplot(2,1,1)
# qwn = np.sqrt(test.qxn**2+test.qyn**2)
# plt.imshow(qwn)

# do iterations
plt.figure()
start_inds, new_inds = test.run_water_iteration()
for j in range(0,len(start_inds)):
    plt.scatter(start_inds[j][1],start_inds[j][0],c='b',s=7)
    plt.scatter(new_inds[j][1],new_inds[j][0],c='r',s=4)
plt.imshow(params.stage)
plt.show()

# # make post-iteration plot
# plt.subplot(2,1,2)
# qwn = np.sqrt(test.qxn**2+test.qyn**2)
# plt.imshow(qwn)
# plt.show()
