"""Example of the workflow with gridded anuga output data."""

import numpy as np
import os.path
import dorado.particle_track as pt
from dorado.routines import get_state
from dorado.routines import plot_state
import matplotlib.pyplot as plt

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
params.stage = np.copy(depth)  # using depth as standin for stage parameter
params.qx = qx
params.qy = qy
# smaller cell size than the other examples to tighten the range of
# travel time values that are obtained (aka increase travel time resolution)
params.dx = 10.
params.theta = 1.0
params.model = 'Anuga'

# particle info
seed_xloc = list(range(20, 30))
seed_yloc = list(range(48, 53))
Np_tracer = 50

# Apply the parameters to run the model
particle = pt.Particles(params)
np.random.seed(0)
# generate particles
particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
# run model until all particles have travelled for about 1.5 hours
walk_data = particle.run_iteration(target_time=2100)
_, _, finaltimes = get_state(walk_data, iteration=-1)

plot_state(params.depth, walk_data, iteration=0, c='b')
plot_state(params.depth, walk_data, iteration=-1, c='r')

# print target travel time and list of the particle travel times
print('Prescribed target travel time: 2100 seconds')
print('List of particle travel times for final particle locations: ' +
      str(np.round(finaltimes)))

# make plot of initial and final particle locations
plt.title('Initial and Final Particle Locations')
plt.show()
