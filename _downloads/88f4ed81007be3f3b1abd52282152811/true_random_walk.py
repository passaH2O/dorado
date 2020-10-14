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
from dorado import particle_track
from dorado import routines

# Here we define the parameters for the particle routing. Since we are
# simulating the generic random walk here, we set our water depth and stage to
# all 1s, and the flow discharge components to all 0s. In this case we are
# seeding 50 particles in the center of our domain. The weighting parameters
# theta and gamma are set to 0.

# create params and then assign the parameters
params = particle_track.modelParams()

# define the params variables
params.depth = np.ones((100, 100))
params.stage = np.ones((100, 100))
params.qx = np.zeros_like(params.depth)
params.qy = np.zeros_like(params.depth)
params.dx = 50.
params.theta = 0.0
params.gamma = 0.0
params.model = 'None'

# particle info
seed_xloc = list(range(45, 56))
seed_yloc = list(range(45, 56))
Np_tracer = 50

np.random.seed(0)  # fix random seed for the example

# initialize the particles class and generate the initial particles
particles = particle_track.Particles(params)
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

# Now we can apply the high-level API function, `steady_plots`, to simulate the
# movement of the particles we have defined in `params` for 50 iterations.

# using steady (time-invariant) plotting routine
walk_data = routines.steady_plots(particles, 50, 'true_random_walk')
