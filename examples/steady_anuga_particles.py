"""Make an example of the workflow with gridded anuga output data."""

import numpy as np
from particlerouting.routines import steady_plots
from particlerouting.routines import draw_travel_path

# Define the parameters that are being used


# define an empty class
class pobj():
    """Empty class for parameters."""

    pass


# create params and then assign the parameters
params = pobj()

# load some variables from a anuga output so stage is varied
data = np.load('ex_anuga_data.npz')

# pull depth and stage from that data
depth = data['depth']
qx = data['qx']
qy = data['qy']

# define the params variables
params.stage = depth
params.depth = depth
params.qx = qx
params.qy = qy

params.seed_xloc = list(range(20, 30))
params.seed_yloc = list(range(48, 53))
params.Np_tracer = 50
params.dx = 50.
params.theta = 1.0
params.model = 'Anuga'

# Apply the parameters to run the particle routing model
np.random.seed(0)  # fix the random seed
# using steady (time-invariant) plotting routine
walk_data = steady_plots(params, 50, 'steady_anuga_example')

# let's visualize a few of these particle travel paths
draw_travel_path(depth, walk_data, [0,1,2,3],
                 'steady_anuga_example/figs/travel_paths.png')
