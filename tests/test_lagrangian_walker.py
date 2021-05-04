from __future__ import division, print_function, absolute_import
from builtins import range, map
from math import floor, sqrt, pi
import pytest

import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np
import dorado.lagrangian_walker as lw
import dorado.particle_track as pt

# defining some of the geometric/walk values up top
jwalk = np.array([[-1, -1, -1],
                  [0, 0, 0],
                  [1, 1, 1]])

iwalk = np.array([[-1, 0, 1],
                  [-1, 0, 1],
                  [-1, 0, 1]])

distances = np.array([[np.sqrt(2), 1, np.sqrt(2)],
                     [1, 1, 1],
                     [np.sqrt(2), 1, np.sqrt(2)]])

ivec = np.array([[-np.sqrt(0.5), 0, np.sqrt(0.5)],
                 [-1, 0, 1],
                 [-np.sqrt(0.5), 0, np.sqrt(0.5)]])

jvec = np.array([[-np.sqrt(0.5), -1, -np.sqrt(0.5)],
                 [0, 0, 0],
                 [np.sqrt(0.5), 1, np.sqrt(0.5)]])

angles = np.array([[3*pi/4, pi/2, pi/4],
                   [pi, 0, 0],
                   [5*pi/4, 3*pi/2, 7*pi/4]])

# set up a set of model parameters that the test functions can use
params = pt.modelParams()
# define a bunch of expected values
params.stage = np.ones((5, 5))
params.cell_type = np.zeros_like(params.stage)
params.qy = params.stage.copy()
params.qx = np.zeros((5, 5))
params.ivec = ivec
params.jvec = jvec
params.distances = distances
params.dry_depth = 0.1
params.gamma = 0.02
params.theta = 1
params.steepest_descent = True
params.depth = params.stage.copy()
params.depth[2, 2] = 10.0  # define index 8 as the deepest
params.seed_xloc = [1]
params.seed_yloc = [1]
params.Np_tracer = 1
params.dx = 1
params.model = 'DeltaRCM'


# defining the unit tests, one per function in the particle_tools.py
def test_random_pick_seed():
    '''
    Test for function random_pick_seed
    '''
    choices = [0]
    probs = 1
    # should return the only option from the choices input
    assert lw.random_pick_seed(choices, probs) == choices[0]


def test_get_weight():
    '''
    Test for function get_weight
    '''
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # define particles class
    params.depth = params.stage.copy()
    particles = pt.Particles(params)
    # make assertion
    assert lw.get_weight(particles, ind) == 5


def test_calculate_new_ind():
    '''
    Test for function calculate_new_ind within Tools class
    '''
    # assign old index
    old_ind = [1, 1]
    # assign new cell location
    new_cell = 0
    # expect new cell to be in location (0,0)
    assert lw.calculate_new_ind(old_ind, new_cell, iwalk, jwalk) == (0,0)


def test_step_update_straight():
    '''
    Test for function step_update within Tools class
    '''
    # set cell size to 1
    dx = 1.
    # define new cell location
    new_cell = 1
    # expect distance between new and old locations to be 1
    # would expect sqrt(2) if the step was diagonal instead of vertical
    assert lw.step_update(new_cell, distances, dx) == 1


def test_step_update_diagonal():
    '''
    Test for function step_update within Tools class
    '''
    # set cell size to 1
    dx = 1.
    # define new cell location
    new_cell = 2
    # expect distance between new and old locations to be sqrt(2)
    # would expect 1 if the step was vertical instead of diagonal
    assert lw.step_update(new_cell, distances, dx) == sqrt(2)


def test_calc_travel_times():
    '''
    Test for function calc_travel_times within Tools class
    '''
    # define particles
    params.diff_coeff = 0
    params.depth = np.ones_like(params.stage)
    params.qx = np.zeros_like(params.stage)
    params.qy = np.ones_like(params.stage)
    particle = pt.Particles(params)
    # define old ind
    old_ind = [1,1]
    # define new ind
    new_ind = [0,1]
    # hardset angles/velocities to make computation more obvious
    particle.velocity = np.ones_like(particle.velocity)
    particle.velocity[0, 0:2] = 3
    particle.velocity_angle = np.ones_like(particle.velocity_angle)
    # get time
    trav_time = lw.calc_travel_times(particle, 1, old_ind, new_ind, 1)
    # expect to return the value 0.5 (inverse of the avg velocity 2)
    assert trav_time == pytest.approx(0.5609806565385976)


def test_check_for_boundary():
    '''
    Test for function check_for_boundary within Tools class
    '''
    # define padded cell types for tools class
    cell_type = np.ones((3,3))
    # define an edge (type==-1)
    cell_type[0,0:2] = -1
    # define new ind
    new_ind = [[0,1]]
    # define current ind
    current_ind = [[1,1]]
    # expect to return the current ind because proposed new ind is an edge
    assert lw.check_for_boundary(new_ind,current_ind,cell_type) == current_ind


def test_random_pick():
    '''
    Test for function random_pick within Tools class
    '''
    # define probs array of zeros with a single 1 value
    probs = np.zeros((8,))
    probs[0] = 1
    # should return first index
    assert lw.random_pick(probs) == 0


def test_get_weight_norm():
    '''
    Test for function get_weight within Tools class
    '''
    # define particles
    params.qy = np.ones_like(params.stage)
    params.qx = np.zeros_like(params.stage)
    params.stage[1,1] = 100.0
    particles = pt.Particles(params)
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert lw.get_weight(particles, ind) == 5


def test_get_weight_deep():
    '''
    Test for function get_weight within Tools class
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.ones((5,5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5,5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances*np.nan
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.depth[2,2] = 10.0 # define index 8 as the deepest
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert lw.get_weight(particles, ind) == 8


def test_make_weight_deep():
    '''
    Test for function make_weight within lagrangian_walker
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.ones((5, 5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5, 5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances*np.nan
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.depth[2,2] = 10.0  # define index 8 as the deepest
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert lw.get_weight(particles, ind) == 8


def test_make_weight_shallow():
    '''
    Test calculating weights when all neighboring cells are dry
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.zeros((5, 5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5, 5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances*np.nan
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.depth[1, 1] = 10.0  # particle location will be wet
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # make assertions about weights
    # at index, index[4] (self) will be 1 while neighbors will be 0
    assert particles.weight[1, 1, 4] == 1.0
    assert np.sum(particles.weight[1, 1, :]) == 1.0
    # weights at boundary cells should be 0
    assert np.all(np.sum(particles.weight[0, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[-1, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, 0, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, -1, 4]) == 0.0)


def test_make_weight_equal_opportunity():
    '''
    Test calculating weights
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.zeros((5, 5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5, 5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances*np.nan
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.depth[1, 1] = 10.0  # same depth at some neighboring cells
    tools.depth[1, 2] = 10.0  # same depth at some neighboring cells
    tools.depth[2, 2] = 10.0  # same depth at some neighboring cells
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # make assertions about weights
    # at wet locations, 2 neighbors will be equiprobable
    assert np.sum(particles.weight[1, 1, :]) == 3.0
    assert np.sum(particles.weight[1, 2, :]) == 3.0
    assert np.sum(particles.weight[2, 2, :]) == 3.0
    # weights at boundary cells should be 0
    assert np.all(np.sum(particles.weight[0, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[-1, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, 0, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, -1, 4]) == 0.0)


def test_make_weight_unequal_opportunity():
    '''
    Test calculating weights w/ different depths in cells
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.zeros((5, 5))
    tools.cell_type = tools.stage.copy()
    tools.qy = tools.stage.copy()
    tools.qx = tools.stage.copy()
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.depth[1, 1] = 10.0  # same depth at some neighboring cells
    tools.depth[1, 2] = 5.0  # less depth at some neighboring cells
    tools.depth[2, 2] = 5.0  # less depth at some neighboring cells
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # make assertions about weights
    # at index, staying put index[4] higher probability than neighbors
    assert particles.weight[1, 1, 4] > particles.weight[1, 1, 5]
    assert particles.weight[1, 1, 4] > particles.weight[1, 1, 8]
    # but the two neighbors should be equiprobable
    assert particles.weight[1, 1, 5] == particles.weight[1, 1, 8]
    # weights at boundary cells should be 0
    assert np.all(np.sum(particles.weight[0, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[-1, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, 0, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, -1, 4]) == 0.0)


def test_wet_boundary_no_weight():
    '''
    Confirm that even if the boundary cell is deep enough, its weight == 0.
    '''
    tools = pt.modelParams()
    # define a bunch of expected values
    tools.stage = np.ones((5, 5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5, 5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances*np.nan
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy() * 10.0
    tools.seed_xloc = [1]
    tools.seed_yloc = [1]
    tools.Np_tracer = 1
    tools.dx = 1
    # define particles
    particles = pt.Particles(tools)
    # assert weights at boundary cells should be 0
    assert np.all(np.sum(particles.weight[0, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[-1, :, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, 0, 4]) == 0.0)
    assert np.all(np.sum(particles.weight[:, -1, 4]) == 0.0)
    # assert that weights everywhere else are not 0
    assert np.all(np.sum(particles.weight[1:-1, 1:-1, 4]) != 0.0)
    # assert that depths everywhere are 10.0
    assert np.all(particles.depth == 10.0)
