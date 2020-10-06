"""Make an example of the workflow with deltarcm output data."""

import numpy as np
import os.path
from dorado.routines import steady_plots
import dorado.particle_track as pt

# load some variables from a deltarcm output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_deltarcm_data.npz')
data = np.load(data_path)

# pull depth and stage from that data
stage = data['stage']
depth = data['depth']

# create params and then assign the parameters
params = pt.modelParams()

# define the params variables
params.stage = stage
params.depth = depth
params.dx = 50.

# don't have any discharge data so we use zeros in place of it
# weighting scheme uses both discharge and depth so can still weight
# based on the water depths
params.qx = np.zeros(np.shape(params.depth))
params.qy = np.zeros(np.shape(params.depth))
params.theta = 1.0
params.model = 'DeltaRCM'  # say that our inputs are from DeltaRCM

# define initial particle locations and number of particles
seed_xloc = list(range(15, 17))
seed_yloc = list(range(137, 140))
Np_tracer = 50

# fix the random seed for the example
np.random.seed(0)

# initialize the particles object and generate some particles
particles = pt.Particles(params)
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

# using steady (time-invariant) plotting routine
walk_data = steady_plots(particles, 50, 'steady_deltarcm_example')
