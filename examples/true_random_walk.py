# example of turning off the weights and reducing problem to a true random walk

import numpy as np
from particlerouting.routines import steady_plots

### Define the parameters that are being used

# define an empty class
class pobj():
    pass

# create params and then assign the parameters
params = pobj()

# define the params variables
params.depth = np.ones((100,100))
params.qx = np.zeros_like(params.depth)
params.qy = np.zeros_like(params.depth)

params.seed_xloc = list(range(45,56))
params.seed_yloc = list(range(45,56))
params.Np_tracer = 50
params.dx = 50.
params.theta = 0.0
params.gamma = 0.0
params.model = 'None'

### Apply the parameters to run the particle routing model

# using steady (time-invariant) plotting routine
steady_plots(params, 50, 'true_random_walk')
