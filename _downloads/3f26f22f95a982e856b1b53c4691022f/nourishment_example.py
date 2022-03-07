"""example of the nourishment functions on DeltaRCM topography"""
import numpy as np
import matplotlib.pyplot as plt
import dorado
import os.path

# load some variables from a deltarcm output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_deltarcm_data.npz')
data = np.load(data_path)

# create params and then assign the parameters
params = dorado.particle_track.modelParams()

# define the params variables
params.stage = data['stage']
params.depth = data['depth']
params.qx = data['qx']
params.qy = data['qy']
params.dx = 50.
params.gamma = 0.5
params.model = 'DeltaRCM'  # say that our inputs are from DeltaRCM

# define initial particle locations and number of particles
seed_xloc = list(range(15, 17))
seed_yloc = list(range(137, 140))
Np_tracer = 300

# fix the random seed for the example
np.random.seed(0)

# initialize the particles object and generate some particles
particles = dorado.particle_track.Particles(params)
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

# using steady (time-invariant) plotting routine
walk_data = dorado.routines.steady_plots(particles, 120,
                                         'nourishment_example')

# let's compute the "nourishment area" for this seed location,
# which shows the fraction of particles that visited each cell
# in the domain. this function includes a bit of gaussian
# filtering to smooth out some of the stochasticity
visit_freq = dorado.particle_track.nourishment_area(walk_data,
                                                    params.depth.shape)

# visualize nourishment area using built-in plotting routine
ax = dorado.routines.show_nourishment_area(visit_freq,
                                           params.depth, walk_data)
plt.savefig(os.path.join('nourishment_example', 'NourishmentArea.png'),
            bbox_inches='tight')

# now, let's look at the "nourishment time," i.e. how long on
# average each particle spends in a given cell once it gets there
mean_times = dorado.particle_track.nourishment_time(walk_data,
                                                    params.depth.shape,
                                                    clip=95)

# and again use the built-in plotting routine
ax = dorado.routines.show_nourishment_time(mean_times, params.depth,
                                           walk_data)
plt.savefig(os.path.join('nourishment_example', 'NourishmentTime.png'),
            bbox_inches='tight')
