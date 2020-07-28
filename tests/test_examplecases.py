from __future__ import division, print_function, absolute_import
from builtins import range, map
from math import floor, sqrt, pi
import pytest
import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np
import dorado.particle_track as pt
from dorado.particle_track import Particle

# set up DeltaRCM test case
rcm_data = np.load('tests/data/ex_deltarcm_data.npz')
params = pt.params()
params.stage = rcm_data['stage']
params.depth = rcm_data['depth']
# just do one particle in a set location for this test
params.seed_xloc = [16]
params.seed_yloc = [140]
params.Np_tracer = 1
params.dx = 50.
# no discharge data so use 0s as substitute
params.qx = np.zeros_like(params.depth)
params.qy = np.zeros_like(params.depth)
params.model = 'DeltaRCM'


def test_few_steps_RCM():
    '''
    Test running a few steps
    '''
    # defining / initializing
    particle = Particle(params) # define the particle
    np.random.seed(0) # fix the random seed
    all_walk_data = None # init the walk data

    # 3 iterations
    for i in list(range(0,3)):
        all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data)

    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # check all positions and times
    assert pytest.approx(all_walk_data['xinds'][0] == [16, 17, 18, 19])
    assert pytest.approx(all_walk_data['yinds'][0] == [140, 141, 141, 141])
    assert pytest.approx(all_walk_data['travel_times'][0] == [0.0,
                                                              1414213.56594,
                                                              2414213.56846,
                                                              3414213.57099])


def test_set_time_RCM_previousdata():
    '''
    Test setting a time target when using old walk data
    '''
    # defining / initializing
    particle = Particle(params) # define the particle
    np.random.seed(0) # fix the random seed

    # iterate once to generate previous walk data
    all_walk_data = particle.run_iteration()

    # set time
    all_walk_data = particle.run_iteration(previous_walk_data=all_walk_data,
                                           target_time=5e6)
    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # check all positions and times
    assert pytest.approx(all_walk_data['xinds'][0] == [16, 17, 18, 19, 20])
    assert pytest.approx(all_walk_data['yinds'][0] == [140, 141, 141, 141, 142])
    assert pytest.approx(all_walk_data['travel_times'][0] == [0.0,
                                                              1414213.56594,
                                                              2414213.56846,
                                                              3414213.57099,
                                                              4828427.13693])


def test_set_time_RCM():
    '''
    Test setting a time target
    '''
    # defining / initializing
    particle = Particle(params) # define the particle
    np.random.seed(0) # fix the random seed

    # set time
    all_walk_data = particle.run_iteration(target_time=5e6)

    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # check all positions and times
    assert pytest.approx(all_walk_data['xinds'][0] == [16, 17, 18, 19, 20])
    assert pytest.approx(all_walk_data['yinds'][0] == [140, 141, 141, 141, 142])
    assert pytest.approx(all_walk_data['travel_times'][0] == [0.0,
                                                              1414213.56594,
                                                              2414213.56846,
                                                              3414213.57099,
                                                              4828427.13693])

# set up anuga test case - going to use subset of domain
an_data = np.load('tests/data/ex_anuga_data.npz')
an_params = pt.params()
an_params.stage = an_data['depth'][40:61,40:61] # just using depth here
an_params.depth = an_data['depth'][40:61,40:61]
# just do one particle in a set location for this test
an_params.seed_xloc = [5]
an_params.seed_yloc = [10]
an_params.Np_tracer = 1
an_params.dx = 50.
# no discharge data so use 0s as substitute
an_params.qx = an_data['qx'][40:61,40:61]
an_params.qy = an_data['qy'][40:61,40:61]
an_params.model = 'Anuga'


def test_few_steps_anuga():
    '''
    Test running a few steps
    '''
    # defining / initializing
    an_particle = Particle(an_params) # define the particle
    np.random.seed(0) # fix the random seed
    all_walk_data = None # init the walk data

    # 3 iterations
    for i in list(range(0,3)):
        all_walk_data = an_particle.run_iteration(previous_walk_data=all_walk_data)

    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == an_params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == an_params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # check all positions and times
    assert pytest.approx(all_walk_data['xinds'][0] == [5, 6, 7, 8])
    assert pytest.approx(all_walk_data['yinds'][0] == [10, 10, 10, 10])
    assert pytest.approx(all_walk_data['travel_times'][0] == [0.0,
                                                              2.68302,
                                                              5.27901,
                                                              7.83891])

def test_boundary_anuga():
    '''
    Test running into the boundary
    '''
    # defining / initializing
    an_particle = Particle(an_params) # define the particle
    np.random.seed(0) # fix the random seed
    all_walk_data = None # init the walk data

    # 20 iterations
    for i in list(range(0,20)):
        all_walk_data = an_particle.run_iteration(previous_walk_data=all_walk_data)

    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == an_params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == an_params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # particle reaches boundary after 18 iterations
    # data shouldn't be recorded after boundary is reached
    assert len(all_walk_data['xinds'][0]) == 15
    assert len(all_walk_data['yinds'][0]) == 15
    assert len(all_walk_data['travel_times'][0]) == 15


def test_boundary_travel_time_anuga():
    '''
    Test running into the boundary and not reaching travel time target
    '''
    # defining / initializing
    an_particle = Particle(an_params) # define the particle
    np.random.seed(0) # fix the random seed

    # set target time for iterations
    all_walk_data = an_particle.run_iteration(target_time=1000.0)

    # make assertions
    # check initial position and travel time
    assert all_walk_data['xinds'][0][0] == an_params.seed_xloc[0]
    assert all_walk_data['yinds'][0][0] == an_params.seed_yloc[0]
    assert all_walk_data['travel_times'][0][0] == 0.0
    # particle reaches boundary after 18 iterations
    # data shouldn't be recorded after boundary is reached
    assert len(all_walk_data['xinds'][0]) == 9
    assert len(all_walk_data['yinds'][0]) == 9
    assert len(all_walk_data['travel_times'][0]) == 9
