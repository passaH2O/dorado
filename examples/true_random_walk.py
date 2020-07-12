"""
True Random Walk
================

In this example the routing weights are turned off and the problem is reduced
to a true random walk.

"""

# First we will import our libraries, numpy and matplotlib are used to handle
# assigning particle parameters and plotting the results.
# The particle_track module is loaded to help define the particle parameters.
# The routines module is used to access the high-level API functions for
# particle movement on a steady flow field, and plotting particle locations.

import numpy as np
import matplotlib.pyplot as plt
from particlerouting import particle_track
from particlerouting import routines

# Here we define the parameters for the particle routing. Since we are
# simulating the generic random walk here, we set our water depth and stage to
# all 1s, and the flow discharge components to all 0s. In this case we are
# seeding 50 particles in the center of our domain. The weighting parameters
# theta and gamma are set to 0.

# create params and then assign the parameters
params = particle_track.params()

# define the params variables
params.depth = np.ones((100, 100))
params.stage = np.ones((100, 100))
params.qx = np.zeros_like(params.depth)
params.qy = np.zeros_like(params.depth)

params.seed_xloc = list(range(45, 56))
params.seed_yloc = list(range(45, 56))
params.Np_tracer = 50
params.dx = 50.
params.theta = 0.0
params.gamma = 0.0
params.model = 'None'

# Now we can apply the high-level API function, `steady_plots`, to simulate the
# movement of the particles we have defined in `params` for 50 iterations.
# Default behavior for this function is to write and save figures and particle
# location data to the disk, for this example we are turning that option off.

# using steady (time-invariant) plotting routine
walk_data = routines.steady_plots(params, 50, 'true_random_walk',
                                  save_output=False)
