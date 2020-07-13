"""Plot for quickstart demo 2 particle locations."""
import matplotlib.pyplot as plt
from particlerouting.example_data.define_params import make_anuga_params
import particlerouting as pr
import numpy as np

np.random.seed(1)  # make result consistent for docs
anugaparams = make_anuga_params()
particle = pr.particle_track.Particle(anugaparams)
walk_data = particle.run_iteration(target_time=2100)
plt.figure(figsize=(8, 5))
pr.routines.plot_initial(anugaparams.depth, walk_data)
pr.routines.plot_final(anugaparams.depth, walk_data)
plt.title('Initial and Final Particle Locations')
plt.tight_layout()
plt.show()
