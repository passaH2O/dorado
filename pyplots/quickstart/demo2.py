"""Plot for quickstart demo 2 particle locations."""
import matplotlib.pyplot as plt
from dorado.example_data.define_params import make_anuga_params
import dorado as pr
import numpy as np

np.random.seed(1)  # make result consistent for docs
anugaparams = make_anuga_params()
particle = pr.particle_track.Particle(anugaparams)
walk_data = particle.run_iteration(target_time=2100)
plt.figure(figsize=(8, 5))
pr.routines.plot_state(anugaparams.depth, walk_data, iteration=0, c='b')
pr.routines.plot_state(anugaparams.depth, walk_data, iteration=-1, c='r')
plt.title('Initial and Final Particle Locations')
plt.tight_layout()
plt.show()
