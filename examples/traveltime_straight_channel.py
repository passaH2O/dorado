"""Example case for particle travel times in a straight channel."""
import numpy as np
import matplotlib.pyplot as plt
import dorado.particle_track as pt

# fix the random seed so it stays the same as weights change
np.random.seed(1)

# create synthetic domain and flow field
domain = np.zeros((100, 50))
depth = np.zeros_like(domain)
stage = np.zeros_like(domain)
u = np.zeros_like(domain)
v = np.zeros_like(domain)
dx = 50.
Np_tracer = 500
seed_xloc = [10]
seed_yloc = [25]

# set up straight channel
depth[:, 10:40] = 1.0
stage[:, 10:40] = 1.0
v[:, 10:40] = -10.0

# choose number of iterations for particle to route
num_iter = 100

# define your 'known' or 'expected' travel time for this simple geometry
# picking expected time from location x=10 to x=70
# (really the boundary of row 70, so 1/2 a cell)
# 59.5 cells * 50 m/cell / 10 m/s = 297.5 seconds
target_row = 70
expected_time = 297.5

# assign particle parameters
params = pt.modelParams()
params.depth = depth
params.stage = stage
params.u = u
params.v = v
params.dx = dx

# set-up figure
plt.figure()
plt.imshow(np.sqrt(u**2 + v**2))
plt.colorbar()
plt.scatter(seed_yloc, seed_xloc, c='k', marker='o', s=5)
# plot the target line where time is measured
plt.plot(np.linspace(0, 50, 100), np.ones(100)*target_row, c='red')
plt.title('Velocity Field')
plt.legend(labels=['Target Row to Measure Times',
                   'Particle Seeding Location'],
           loc='best')
plt.tight_layout()
plt.show()

# do the routing twice, once without any diffusivity added to the travel times
# (diff_coeff==0) then a second time with significant diffusion (diff_coeff==1)
for dc in list(range(0, 2)):
    # set diff_coeff
    if dc == 0:
        params.diff_coeff = 0.0
    else:
        params.diff_coeff = 1.0

    # make particle
    particle = pt.Particles(params)

    # walk it
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    for i in list(range(0, num_iter)):
        walk_data = particle.run_iteration()

    # get travel times associated with particles when they are at coord x=70
    # use the exposure_time function to measure this
    roi = np.zeros_like(depth, dtype='int')
    roi[0:target_row, :] = 1
    target_times = pt.exposure_time(walk_data, roi)

    # plot histogram
    plt.subplot(1, 2, dc+1)
    n, bins, _ = plt.hist(target_times, bins=100, range=(200, 400),
                          histtype='bar', density=True,
                          color=[0.5, 0.5, 1, 0.5])

    # plot expected travel time to row 70
    plt.scatter(expected_time, np.max(n),
                s=75, c='green', marker='x', linewidths=20)

    plt.legend(['Expected Travel Time',
                'Histogram of Final Travel Times'], ncol=2,
               loc='upper left', bbox_to_anchor=(0.0, -0.06), fontsize=16)

    plt.title('Travel Time Distribution at Target Row \n'
              'Diffusion Coefficient : ' + str(params.diff_coeff), fontsize=20)
    plt.xlabel('Travel Time at Target Row [s]', fontsize=16)
    plt.ylabel('Probability Density', fontsize=16)

plt.show()
