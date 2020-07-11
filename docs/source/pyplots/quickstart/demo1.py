"""Plot for quickstart demo 1."""
import matplotlib.pyplot as plt
from particlerouting import routines
from particlerouting.example_data.define_params import make_rcm_params
import numpy as np

np.random.seed(1)  # fix random seed so example is always the same
params = make_rcm_params()
all_walk_data = routines.steady_plots(params, 50, 'demo-1', save_imgs=False)

plt.figure(figsize=(8, 5))
for k in list(range(0, params.Np_tracer)):
    plt.scatter(all_walk_data['yinds'][k][0],
                all_walk_data['xinds'][k][0],
                c='b',
                s=4.0)
    plt.scatter(all_walk_data['yinds'][k][-1],
                all_walk_data['xinds'][k][-1],
                c='r',
                s=4.0)
plt.imshow(params.depth)
plt.title('Depth - Particle Iteration 49')
cbar = plt.colorbar(fraction=0.025, pad=0.03)
cbar.set_label('Water Depth [m]', labelpad=10.0)
plt.axis('scaled')
plt.show()
