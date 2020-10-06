"""Unsteady flow example using gridded anuga output data."""
# Note: This example must be run from within the 'examples' directory!!!

from dorado.routines import unsteady_plots

# set info not contained in saved data
dx = 5.
Np_tracer = 50
seed_xloc = list(range(5, 16))
seed_yloc = list(range(48, 53))

# then apply the unsteady_plots function and let it assign domain parameters
walk_data = unsteady_plots(dx, Np_tracer, seed_xloc, seed_yloc,
                           26, 75., 'unsteady_data',
                           'csv', 'unsteady_output')
