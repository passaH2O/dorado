import pytest

import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np
from particlerouting.particle_tools import Tools

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



# defining the unit tests, one per function in the particle_tools.py
def test_random_pick_seed():
    '''
    Test for function random_pick_seed within Tools class
    '''
    choices = [0]
    probs = 1
    tools = Tools()
    # should return the only option from the choices input
    assert tools.random_pick_seed(choices,probs) == choices[0]



def test_get_weight():
    '''
    Test for function get_weight within Tools class
    '''
    tools = Tools()
    # define a bunch of expected values
    tools.stage = np.ones((5,5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.stage = tools.stage.copy()
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5,5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert tools.get_weight(ind) == 5



def test_calculate_new_ind():
    '''
    Test for function calculate_new_ind within Tools class
    '''
    tools = Tools()
    # assign walk directions
    tools.iwalk = iwalk
    tools.jwalk = jwalk
    # assign old index
    old_ind = [1,1]
    # assign new cell location
    new_cell = 0
    # expect new cell to be in location (0,0)
    assert tools.calculate_new_ind(old_ind,new_cell) == (0,0)



def test_step_update():
    '''
    Test for function step_update within Tools class
    '''
    tools = Tools()
    # define walk directions
    tools.iwalk = iwalk
    tools.jwalk = jwalk
    # define little qxn qyn regions of zeros
    tools.qxn = np.zeros_like(iwalk)
    tools.qyn = np.zeros_like(iwalk)
    # define old index
    old_ind = [1,1]
    # define new index
    new_ind = [0,1]
    # define new cell location
    new_cell = 1
    # expect distance between new and old locations to be 1
    # would expect sqrt(2) if the step was diagonal instead of vertical
    assert tools.step_update(old_ind,new_ind,new_cell) == 1



def test_calc_travel_times():
    '''
    Test for function calc_travel_times within Tools class
    '''
    tools = Tools()
    # define cell size
    tools.dx = 1
    # define some velocities for tools class
    tools.velocity = np.ones((3,3))
    # define faster velocities so the averaging does something
    tools.velocity[0,0:2] = 3
    # define old ind
    old_ind = [1,1]
    # define new ind
    new_ind = [0,1]
    # expect to return the value 0.5 (inverse of the avg velocity 2)
    assert tools.calc_travel_times(old_ind,new_ind,1) == 0.5



def test_check_for_boundary():
    '''
    Test for function check_for_boundary within Tools class
    '''
    tools = Tools()
    # define padded cell types for tools class
    tools.cell_type = np.ones((3,3))
    # define an edge (type==-1)
    tools.cell_type[0,0:2] = -1
    # define new ind
    new_ind = [[0,1]]
    # define current ind
    current_ind = [[1,1]]
    # expect to return the current ind because proposed new ind is an edge
    assert tools.check_for_boundary(new_ind,current_ind) == current_ind



def test_random_pick():
    '''
    Test for function random_pick within Tools class
    '''
    tools = Tools()
    # define probs array of zeros with a single 1 value
    probs = np.zeros((8,))
    probs[0] = 1
    # should return first index
    assert tools.random_pick(probs) == 0



def test_get_weight_norm():
    '''
    Test for function get_weight within Tools class
    '''
    tools = Tools()
    # define a bunch of expected values
    tools.stage = np.ones((5,5))
    tools.cell_type = np.zeros_like(tools.stage)
    tools.qy = tools.stage.copy()
    tools.qx = np.zeros((5,5))
    tools.ivec = ivec
    tools.jvec = jvec
    tools.distances = distances
    tools.dry_depth = 0.1
    tools.gamma = 0.02
    tools.theta = 1
    tools.steepest_descent = True
    tools.depth = tools.stage.copy()
    tools.stage[1,1] = 100.0 # make stage at ind large so it is normalized
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert tools.get_weight(ind) == 5



def test_get_weight_deep():
    '''
    Test for function get_weight within Tools class
    '''
    tools = Tools()
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
    # set the current index
    ind = (1,1)
    # set seed
    np.random.seed(0)
    # make assertion
    assert tools.get_weight(ind) == 8
