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

params.seed_xloc = [50,51,52]
params.seed_yloc = [50,51,52]
params.Np_tracer = 3
params.dx = 50.
params.depth = np.ones((100,100))
params.stage = np.ones((100,100))
params.qx = np.zeros((100,100))
params.qy = np.zeros((100,100))
params.theta = 1.0

# create some discharge
params.qx[40:60,40:60] = 0.25
params.qy[30:45,30:45] = 0.50

# try running it
from particle_track import Particle

test = Particle(params)

test.init_water_iteration()

# make pre-iteration plot
plt.figure()
plt.subplot(2,1,1)
qwn = np.sqrt(test.qxn**2+test.qyn**2)
plt.imshow(qwn)

# do iteration
test.run_water_iteration()

# make post-iteration plot
plt.subplot(2,1,2)
qwn = np.sqrt(test.qxn**2+test.qyn**2)
plt.imshow(qwn)
plt.show()
