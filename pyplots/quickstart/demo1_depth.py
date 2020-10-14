"""Plot for quickstart demo 1."""
import matplotlib.pyplot as plt
from dorado.example_data.define_params import make_rcm_particles

particles = make_rcm_particles()
fig = plt.figure()
ax = plt.gca()
im = ax.imshow(particles.depth)
plt.title('Water Depth')
cax = fig.add_axes([ax.get_position().x1+0.01,
                    ax.get_position().y0,
                    0.02,
                    ax.get_position().height])
cbar = plt.colorbar(im, cax=cax)
plt.show()
