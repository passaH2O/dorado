import pytest

from particlerouting import particle_track
import numpy as np


# init some parameters
class pobj():
    pass

params = pobj
params.seed_xloc = [1]
params.seed_yloc = [1]
params.Np_tracer = 1
params.dx = 1
params.depth = np.ones((3,3))
params.stage = np.ones((3,3))
params.qx = np.zeros((3,3))
params.qy = np.ones((3,3))
params.theta = 1
particle = particle_track.Particle(params)

# testing of the Particle __init__ functionality
def test_xseed():
    assert particle.seed_xloc == params.seed_xloc

def test_yseed():
    assert particle.seed_yloc == params.seed_yloc

def test_Np_tracer():
    assert particle.Np_tracer == params.Np_tracer

def test_dx():
    assert particle.dx == params.dx

def test_depth():
    assert np.all(particle.depth) == np.all(params.depth)

def test_stage():
    assert np.all(particle.stage) == np.all(params.stage)

def test_qx():
    assert np.all(particle.qx) == np.all(params.qx)

def test_qy():
    assert np.all(particle.qy) == np.all(params.qy)

def test_velocity():
    vel = np.sqrt(params.qx**2+params.qy**2)/params.depth
    assert np.all(particle.velocity) == np.all(vel)

def test_theta():
    assert particle.theta == params.theta

def test_itmax():
    assert particle.itmax == 1

def test_dry_depth():
    assert particle.dry_depth == 0.1

def test_gamma():
    assert particle.gamma == 0.05

def test_cell_type():
    assert particle.cell_type[0,0] == 0

def test_distances():
    assert particle.distances[1,1] == 1

def test_ivec():
    assert particle.ivec[1,1] == 0

def test_jvec():
    assert particle.jvec[1,1] == 0

def test_iwalk():
    assert particle.iwalk[1,1] == 0

def test_jwalk():
    assert particle.jwalk[1,1] == 0

def test_qxn():
    assert np.sum(particle.qxn) == 0

def test_qyn():
    assert np.sum(particle.qyn) == 0

def test_sfc_visit():
    assert np.sum(particle.sfc_visit) == 0

def test_sfc_sum():
    assert np.sum(particle.sfc_sum) == 0

def test_pad_stage():
    assert np.all(particle.pad_stage) == np.all(np.pad(params.stage, 1, 'edge'))

def test_pad_depth():
    assert np.all(particle.pad_depth) == np.all(np.pad(params.depth, 1, 'edge'))

def test_pad_cell_type():
    assert np.all(particle.pad_cell_type) == np.all(np.pad(np.zeros_like(params.stage), 1, 'constant', constant_values = -1))

def test_steep():
    assert particle.steepest_descent == False

def test_steep_true():
    params.steepest_descent = True
    particle = particle_track.Particle(params)
    assert particle.steepest_descent == True

def test_steep_other():
    params.steepest_descent = 'other'
    particle = particle_track.Particle(params)
    assert particle.steepest_descent == False



# testing of the run_iteration function
def test_start_pairs_X():
    start_pairs, new_inds, travel_times = particle.run_iteration()
    assert start_pairs[0][0] == params.seed_xloc[0]

def test_start_pairs_Y():
    start_pairs, new_inds, travel_times = particle.run_iteration()
    assert start_pairs[0][1] == params.seed_yloc[0]

def test_start_pairs_X2():
    start_pairs, new_inds, travel_times = particle.run_iteration(start_xindices=params.seed_xloc)
    assert start_pairs[0][0] == params.seed_xloc[0]

def test_start_pairs_Y2():
    start_pairs, new_inds, travel_times = particle.run_iteration(start_yindices=params.seed_yloc)
    assert start_pairs[0][1] == params.seed_yloc[0]

def test_travel_time():
    start_pairs, new_inds, travel_times = particle.run_iteration()
    assert travel_times[0] == 1.0

def test_travel_time_given():
    start_pairs, new_inds, travel_times = particle.run_iteration(start_times=[0.0])
    assert travel_times[0] == 1.0
