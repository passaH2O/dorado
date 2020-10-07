"""Make an example of the workflow with gridded Anuga output data."""

import numpy as np
import os.path
import time
from dorado.parallel_routing import parallel_routing

# for serial run comparison import the regular iterator
from dorado.particle_track import Particles
import dorado.particle_track as pt

# create params and then assign the parameters
params = pt.modelParams()

# load some variables from an anuga output so stage is varied
f_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(f_path, 'ex_anuga_data.npz')
data = np.load(data_path)

# pull depth and stage from that data
depth = data['depth']
qx = data['qx']
qy = data['qy']

# define the params variables
params.depth = depth
params.stage = np.copy(depth)  # use depth as proxy for stage in this example
params.qx = qx
params.qy = qy
params.dx = 50.
params.theta = 1.0
params.model = 'Anuga'

# for parallel routing we define the particle information
seed_xloc = list(range(20, 30))
seed_yloc = list(range(48, 53))
Np_tracer = 200

# Apply the parameters to run the particle routing model

# use 2 cores to route in parallel
print('start parallel')
start_par_time = time.time()
particles = Particles(params)
par_result = parallel_routing(particles, 50, Np_tracer, seed_xloc, seed_yloc, 2)
par_time = time.time() - start_par_time
print('end parallel')

# compare to a serial run
print('start serial')
start_serial_time = time.time()
# do twice to match number of particles parallel is doing
for z in list(range(0, 2)):
    particle = Particles(params)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    # do 50 iterations to match parallel
    for i in list(range(0, 50)):
        all_walk = particle.run_iteration()

# get time
serial_time = time.time() - start_serial_time
print('end serial')

# print times elapsed
print('Serial Compute Time: ' + str(serial_time))
print('Parallel Compute Time: ' + str(par_time))
