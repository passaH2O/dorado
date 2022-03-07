"""Plot for quickstart demo 2 particle locations."""
import matplotlib.pyplot as plt
from dorado.example_data.define_params import make_anuga_params
from dorado.particle_track import Particles
import dorado as pr
import numpy as np

np.random.seed(1)  # make result consistent for docs
params = make_anuga_params()
particles = Particles(params)
seed_xloc = list(range(20, 30))
seed_yloc = list(range(48, 53))
Np_tracer = 50
particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)
walk_data = particles.run_iteration(target_time=2100)
plt.figure(figsize=(8, 5))
pr.routines.plot_state(particles.depth, walk_data, iteration=0, c='b')
pr.routines.plot_state(particles.depth, walk_data, iteration=-1, c='r')
plt.title('Initial and Final Particle Locations')
plt.tight_layout()
plt.show()
