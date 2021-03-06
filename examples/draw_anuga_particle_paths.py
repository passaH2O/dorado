"""Example of method to draw the particle travel paths."""
# Can only be run after steady_anuga_particles.py has been successfully run
import numpy as np
import os
import os.path
import json
from dorado.routines import draw_travel_path

# load the depth data
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_anuga_data.npz')
data = np.load(data_path)
depth = data['depth']

# load the walk data
all_walk_data = json.load(open('steady_anuga_example'+os.sep+'data'+os.sep+'data.txt'))

# Draw the travel path
draw_travel_path(depth, all_walk_data, [0, 1, 2, 3],
                 'steady_anuga_example'+os.sep+'figs'+os.sep+'travel_paths.png',
                 interval=2, plot_legend=True)
