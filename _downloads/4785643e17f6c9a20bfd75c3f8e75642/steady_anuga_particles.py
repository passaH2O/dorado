"""Make an example of the workflow with gridded anuga output data."""

import numpy as np
import os.path
from dorado.routines import steady_plots
import dorado.particle_track as pt

# load some variables from a anuga output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_anuga_data.npz')
data = np.load(data_path)

# pull depth and stage from that data
depth = data['depth']
qx = data['qx']
qy = data['qy']

# create the parameters object and then assign the values
params = pt.modelParams()

# define the params variables
params.stage = np.copy(depth)  # for this example we don't have stage data
params.depth = depth
params.qx = qx
params.qy = qy
params.dx = 50.
params.model = 'Anuga'

# define initial particle locations and number of particles
seed_xloc = list(range(20, 30))
seed_yloc = list(range(48, 53))
Np_tracer = 50

# Apply the parameters to run the particle routing model
np.random.seed(0)  # fix the random seed for example

# initialize particles and generate particles
particles = pt.Particles(params)
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

# using steady (time-invariant) plotting routine
walk_data = steady_plots(particles, 50, 'steady_anuga_example')
