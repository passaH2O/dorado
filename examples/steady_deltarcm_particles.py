"""Make an example of the workflow with deltarcm output data."""

import numpy as np
import os.path
from dorado.routines import steady_plots
from dorado.particle_track import modelParams

# load some variables from a deltarcm output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_deltarcm_data.npz')
data = np.load(data_path)

# pull depth and stage from that data
stage = data['stage']
depth = data['depth']

# create params and then assign the parameters
params = modelParams()

# define the params variables
params.stage = stage
params.depth = depth

params.seed_xloc = list(range(15, 17))
params.seed_yloc = list(range(137, 140))
params.Np_tracer = 50
params.dx = 50.
# don't have any discharge data so we use zeros in place of it
# weighting scheme uses both discharge and depth so can still weight
# based on the water depths
params.qx = np.zeros(np.shape(params.depth))
params.qy = np.zeros(np.shape(params.depth))
params.theta = 1.0
params.model = 'DeltaRCM'  # say that our inputs are from DeltaRCM

# Apply the parameters to run the particle routing model
np.random.seed(0)  # fix the random seed
# using steady (time-invariant) plotting routine
walk_data = steady_plots(params, 50, 'steady_deltarcm_example')
