# Can only be run after steady_deltarcm_particles.py has been successfully run
import numpy as np
from particlerouting.routines import draw_travel_path

### load the depth data
data = np.load('ex_deltarcm_data.npz')
depth = data['depth']

### load the walk data
walkdata = np.load('steady_deltarcm_example/data/data.npz')
all_walk_data = walkdata['all_walk_data']

### Draw the travel path
draw_travel_path(depth, all_walk_data, [0,1,2,3],
                 'steady_deltarcm_example/data/travel_paths.png')
