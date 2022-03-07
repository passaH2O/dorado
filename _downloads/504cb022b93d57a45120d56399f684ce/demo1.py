"""Plot for quickstart demo 1."""
import matplotlib.pyplot as plt
from dorado import routines
from dorado.example_data.define_params import make_rcm_particles
import numpy as np

np.random.seed(1)  # fix random seed so example is always the same
particles = make_rcm_particles()
all_walk_data = routines.steady_plots(particles, 50, 'demo-1',
                                      save_output=False)

fig = plt.figure(figsize=(8, 5))
for k in list(range(0, particles.Np_tracer)):
    plt.scatter(all_walk_data['yinds'][k][0],
                all_walk_data['xinds'][k][0],
                c='b',
                s=4.0)
    plt.scatter(all_walk_data['yinds'][k][-1],
                all_walk_data['xinds'][k][-1],
                c='r',
                s=4.0)
ax = plt.gca()
im = ax.imshow(particles.depth)
plt.title('Depth - Particle Iteration 49')
cax = fig.add_axes([ax.get_position().x1+0.01,
                    ax.get_position().y0,
                    0.02,
                    ax.get_position().height])
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Water Depth [m]', labelpad=10.0)
plt.show()
