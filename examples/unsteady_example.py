"""Unsteady flow example using gridded anuga output data."""

from particlerouting.routines import unsteady_plots
import particlerouting.particle_track as pt

# initialize a parameters object
params = pt.params()

# give params information not contained in the grid
params.dx = 5.
params.Np_tracer = 50
params.seed_xloc = list(range(5, 16))
params.seed_yloc = list(range(48, 53))

# then apply the unsteady_plots function and let it assign domain parameters
unsteady_plots(params, 26, 75., 'unsteady_data',
               'csv', 'unsteady_output')
