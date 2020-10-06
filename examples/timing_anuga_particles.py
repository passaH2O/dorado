"""Example of the workflow with gridded anuga output data"""

import numpy as np
import os.path
from dorado.routines import time_plots
import dorado.particle_track as pt

# load some variables from an anuga output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_anuga_data.npz')
data = np.load(data_path)

# pull depth and stage from that data
depth = data['depth']
qx = data['qx']
qy = data['qy']

# create params and then assign the parameters
params = pt.modelParams()

# define the params variables
params.depth = depth
params.stage = depth  # using depth as stand-in for stage in this example
params.qx = qx
params.qy = qy
params.dx = 50.
params.theta = 1.0
params.model = 'Anuga'

# define particles
seed_xloc = list(range(20, 30))
seed_yloc = list(range(48, 53))
Np_tracer = 50

np.random.seed(0)  # fix random seed for example
particles = pt.Particles(params)
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

# Apply the parameters to run the particle routing model
# using steady (time-invariant) plotting routine
time_plots(params, 50, 'timing_anuga_example')
