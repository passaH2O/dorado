"""Make parameters for the example data."""

import os
import numpy as np
from .. import particle_track as pt


def make_rcm_particles():
    """Function to make the example RCM particles."""

    params = pt.modelParams()
    path = os.path.join(os.path.dirname(__file__), 'ex_deltarcm_data.npz')
    data = np.load(path)
    # pull depth and stage from that data
    stage = data['stage']
    depth = data['depth']
    # define the params variables
    params.stage = stage
    params.depth = depth
    params.dx = 50.
    # don't have any discharge data so we use zeros in place of it
    # weighting scheme uses both discharge and depth so can still weight
    # based on the water depths
    params.qx = np.zeros(np.shape(params.depth))
    params.qy = np.zeros(np.shape(params.depth))
    params.theta = 1.0
    params.model = 'DeltaRCM'  # say that our inputs are from DeltaRCM

    seed_xloc = list(range(15, 17))
    seed_yloc = list(range(137, 140))
    Np_tracer = 50

    particles = pt.Particles(params)
    particles.generate_particles(Np_tracer, seed_xloc, seed_yloc)

    return particles


def make_anuga_params():
    """Function to make the example ANUGA parameters."""

    params = pt.modelParams()
    path = os.path.join(os.path.dirname(__file__), 'ex_anuga_data.npz')
    data = np.load(path)
    # pull depth and stage from that data
    depth = data['depth']
    qx = data['qx']
    qy = data['qy']
    # define the params variables
    params.stage = np.copy(depth)
    params.depth = depth
    params.qx = qx
    params.qy = qy
    params.dx = 10.
    params.theta = 1.0
    params.model = 'Anuga'

    return params
