"""Plot for quickstart demo 2."""
import matplotlib.pyplot as plt
from dorado.example_data.define_params import make_anuga_params
from mpl_toolkits.axes_grid1 import make_axes_locatable

params = make_anuga_params()
plt.figure(figsize=(8, 5))

plt.subplot(1, 2, 1)
ax = plt.gca()
im = ax.imshow(params.qx)
plt.title('x-component of water discharge')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)

plt.subplot(1, 2, 2)
ax = plt.gca()
im = ax.imshow(params.qy)
plt.title('y-component of water discharge')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.show()
