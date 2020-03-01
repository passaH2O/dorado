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

params.inlet = [50,51,52]
params.Np_water = 1
params.Qp_water = 0.01
params.dx = 50.
params.qxn = np.zeros((100,100))
params.qyn = np.zeros((100,100))
params.qwn = np.zeros((100,100))
params.indices = np.array([[0,0]])
params.looped = np.array([0])
params.itmax = 1
params.depth = np.ones((100,100))
params.free_surf_flag = np.array([0])
params.stage = np.ones((100,100))
params.distances = np.array([[sqrt2, 1, sqrt2],
                           [1, 1, 1],
                           [sqrt2, 1, sqrt2]])
params.qx = np.zeros((100,100))
params.qy = np.zeros((100,100))
params.ivec = np.array([[-sqrt05, 0, sqrt05],
                      [-1, 0, 1],
                      [-sqrt05, 0, sqrt05]])
params.jvec = np.array([[-sqrt05, -1, -sqrt05],
                      [0, 0, 0],
                      [sqrt05, 1, sqrt05]])
params.dry_depth = 0.01
params.gamma = 0.02
params.iwalk = np.array([[-1, 0, 1],
                       [-1, 0, 1],
                       [-1, 0, 1]])
params.jwalk = np.array([[-1, -1, -1],
                       [0, 0, 0],
                       [1, 1, 1]])
params.L = 100
params.W = 100
params.cell_type = np.zeros((100,100))
params.L0 = 1
params.CTR = 50
params.sfc_visit = np.zeros((100,100))
params.sfc_sum = np.zeros((100,100))
params.theta_water = 1.0

# try running it
from particle_track import Particle

test = Particle(params)

test.init_water_iteration()

# make pre-iteration plot
plt.figure()
plt.subplot(2,1,1)
plt.imshow(test.qwn)

# do iteration
test.run_water_iteration()

# make post-iteration plot
plt.subplot(2,1,2)
plt.imshow(test.qwn)
plt.show()
