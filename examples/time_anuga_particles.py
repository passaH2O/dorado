# make an example of the workflow with gridded anuga output data

import numpy as np
import particlerouting.particle_track as pt
import matplotlib.pyplot as plt

### Define the parameters that are being used

# define an empty class
class pobj():
    pass

# create params and then assign the parameters
params = pobj()

# load some variables from a deltarcm output so stage is varied
data = np.load('ex_anuga_data.npz')

# pull depth and stage from that data
depth = data['depth']
qx = data['qx']
qy = data['qy']

# define the params variables
params.depth = depth
params.qx = qx
params.qy = qy

params.seed_xloc = list(range(20,30))
params.seed_yloc = list(range(48,53))
params.Np_tracer = 50
params.dx = 50.
params.theta = 1.0
params.model = 'Anuga'

### Apply the parameters to run the particle routing model
particle = pt.Particle(params)
# run model until all particles have travelled for about 1.5 hours
start_inds, end_inds, travel_times = particle.run_iteration(time_step=5400)

# make plot of initial and final particle locations
plt.figure(figsize=(4,4),dpi=200)
for k in range(0,len(start_inds)):
    plt.scatter(start_inds[k][1],start_inds[k][0],c='b',s=0.75)
    plt.scatter(end_inds[k][1],end_inds[k][0],c='r',s=0.75)
plt.imshow(params.depth)
plt.title('Depth')
plt.colorbar()
plt.show()
