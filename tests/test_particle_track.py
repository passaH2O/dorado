from __future__ import division, print_function, absolute_import
from builtins import range, map
from math import floor, sqrt, pi
import copy
import pytest
import io
import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

from dorado import particle_track
import numpy as np

# init some parameters
class pobj():
    pass

params = pobj
seed_xloc = [1]
seed_yloc = [1]
Np_tracer = 1
params.dx = 1
params.depth = np.ones((3,3))
params.stage = np.ones((3,3))
params.qx = np.zeros((3,3))
params.qy = np.ones((3,3))
params.theta = 1
params.model = 'DeltaRCM'
goodparams = copy.deepcopy(params)  # don't corrupt these good parameters
particle = particle_track.Particles(params)

# testing of the Particle __init__ functionality
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

def test_neg_diffcoeff():
    params.diff_coeff = -1.0
    particle = particle_track.Particles(params)
    assert particle.diff_coeff == 0.0

def test_big_diffcoeff():
    params.diff_coeff = 3.0
    particle = particle_track.Particles(params)
    assert particle.diff_coeff == 3.0

def test_dry_depth():
    assert particle.dry_depth == 0.1

def test_gamma():
    assert particle.gamma == 0.05

def test_cell_type():
    assert particle.cell_type[1,1] == 0
    assert np.all(particle.cell_type[0,:] == [-1,-1,-1])
    assert np.all(particle.cell_type[-1,:] == [-1,-1,-1])
    assert np.all(particle.cell_type[:,0] == [-1,-1,-1])
    assert np.all(particle.cell_type[:,-1] == [-1,-1,-1])

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

def test_steep():
    assert particle.steepest_descent == False

def test_steep_true():
    params.diff_coeff = 'invalid'
    params.steepest_descent = True
    particle = particle_track.Particles(params)
    assert particle.steepest_descent == True

def test_steep_other():
    params.steepest_descent = 'other'
    particle = particle_track.Particles(params)
    assert particle.steepest_descent == False

# testing of the run_iteration function
def test_start_pairs_X():
    particle = particle_track.Particles(params)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    all_walk_data = particle.run_iteration()
    assert all_walk_data['xinds'][0][0] == seed_xloc[0]

def test_start_pairs_Y():
    particle = particle_track.Particles(params)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    all_walk_data = particle.run_iteration()
    assert all_walk_data['yinds'][0][0] == seed_yloc[0]

def test_exact_locations():
    num_ps = 2
    particle = particle_track.Particles(params)
    particle.generate_particles(num_ps, seed_xloc, seed_yloc, method='exact')
    all_walk_data = particle.run_iteration()
    assert all_walk_data['xinds'][0] == seed_xloc
    assert all_walk_data['yinds'][0] == seed_yloc

def test_exact_locations_overflow():
    """Test where number of particles exceeds number of seed locations."""
    num_ps = 3
    particle = particle_track.Particles(params)
    particle.generate_particles(num_ps, [0, 1], [0, 1], method='exact')
    all_walk_data = particle.run_iteration()
    assert len(all_walk_data['xinds']) == num_ps
    assert len(all_walk_data['yinds']) == num_ps

def test_no_explicit_generation():
    """Test reading of Np_tracer from self if some walk_data exists."""
    # create some walk data (would exist from a previous run)
    walk_data = {'xinds': [[1], [1]], 'yinds': [[1], [1]],
                 'travel_times': [[0.0], [0.0]]}
    # init particle class and stick walk data in it
    # init defines Np_tracer as 0, but it is 2 in walk_data
    particle = particle_track.Particles(params)
    particle.walk_data = walk_data
    # assert that Np_tracer is 0
    assert particle.Np_tracer == 0
    # run an iteration and see if Np_tracer is corrected
    particle.run_iteration()
    # assert that number of particles has been correctly identified
    assert particle.Np_tracer == 2

def test_travel_time():
    # Particle doesn't travel in the 3x3 space due to the 'sticky' edge
    # conditions so check that travel time is 0 and particle hasn't moved
    particle = particle_track.Particles(params)
    np.random.seed(0)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    all_walk_data = particle.run_iteration()
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0

def test_init_params():
    # test initialization of params class
    params = particle_track.modelParams()
    # make assertions
    assert params.dx is None
    assert params.depth is None
    assert params.stage is None
    assert params.qx is None
    assert params.qy is None
    assert params.u is None
    assert params.v is None


def test_previous_walk_data():
    # test of loading previously defined walk data
    particle = particle_track.Particles(params)
    np.random.seed(0)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    old_walk_data = particle.run_iteration()
    # try to do another walk - test just makes sure code doesn't break
    np.random.seed(0)
    all_walk_data = particle.run_iteration()
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0


def test_generate_twice():
    # test ability to generate particles multiple times
    particle = particle_track.Particles(params)
    np.random.seed(0)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    old_walk_data = particle.run_iteration()
    # try to do another walk - test just makes sure code doesn't break
    np.random.seed(0)
    all_walk_data = particle.run_iteration()
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0
    assert len(all_walk_data['xinds']) == 2


def test_use_walk_data():
    # use walk data without using the generator function
    particle1 = particle_track.Particles(params)
    np.random.seed(0)
    particle1.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    particle1.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    old_walk_data = particle1.run_iteration()
    # create new particles class and try to use that old walk data
    particle2 = particle_track.Particles(params)
    particle2.generate_particles(0, [], [], previous_walk_data=old_walk_data)
    all_walk_data = particle2.run_iteration()
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0
    assert len(all_walk_data['xinds']) == 2


def test_manual_reset():
    # test resetting walk data
    particle = particle_track.Particles(params)
    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
    assert particle.walk_data['xinds'][0][0] == 1
    particle.clear_walk_data()
    assert (particle.walk_data is None) is True


class TestValueErrors:
    """
    Catching the ValueErrors and edge cases
    """

    def test_dx(self):
        params = particle_track.modelParams()
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_broken_depth(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = 'badstring'
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_missing_depth(self):
        params = particle_track.modelParams()
        params.dx = 1
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_depth_via_stage_topo(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.stage = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        particle = particle_track.Particles(params)
        # should work so make assertion
        assert np.all(particle.depth == params.stage-params.topography)

    def test_depth_via_stage_topo_broken(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.stage = 'badstring'
        params.topography = np.zeros((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_stage_broken(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.stage = 'badstring'
        params.depth = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_stage_via_depth_topo(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.topography = np.zeros((3,3))
        params.depth = np.ones((3,3))
        params.depth[0,0] = 0
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        particle = particle_track.Particles(params)
        # should work so make assertion
        assert np.isnan(particle.stage[0,0]) == True
        assert particle.stage[1,1] == 1.

    def test_stage_via_depth_topo_broken(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.topography = 'badstring'
        params.depth = np.ones((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_missing_stage(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = np.ones((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_rcm_model_uv(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        params.model = 'DeltaRCM'
        particle = particle_track.Particles(params)
        # should work so make assertion
        assert np.all(particle.v == params.v)
        assert np.all(particle.u == params.u)
        assert np.all(particle.qx == params.u*params.depth)
        assert np.all(particle.qy == params.v*params.depth)

    def test_rcm_model_uv_broken(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        params.model = 'DeltaRCM'
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_model_q(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.qx = np.ones((3,3))
        params.qy = np.ones((3,3))
        particle = particle_track.Particles(params)
        # should work so make assertion
        assert np.all(particle.v == particle.qy*particle.depth/(particle.depth**2+1e-8))
        assert np.all(particle.u == particle.qx*particle.depth/(particle.depth**2+1e-8))
        assert np.all(particle.qx == -1*params.qy)
        assert np.all(particle.qy == params.qx)

    def test_model_uv_broken(self):
        params = particle_track.modelParams()
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_bad_Np_tracer(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(TypeError):
            particle.generate_particles('invalid', seed_xloc, seed_yloc)

    def test_bad_seed_xloc(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(TypeError):
            particle.generate_particles(1, 'invalid', seed_yloc)

    def test_bad_seed_yloc(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(TypeError):
            particle.generate_particles(1, seed_xloc, 'invalid')

    def test_bad_method(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(ValueError):
            particle.generate_particles(1, seed_xloc, seed_yloc, method='bad')

    def test_invalid_method(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(TypeError):
            particle.generate_particles(1, seed_xloc, seed_yloc, method=5)

    def test_unstruct2grid_break(self):
        coords = [(10.5, 10.1),
                  (10.1, 15.1),
                  (15.2, 20.2)]
        quantity = [1, 2]
        cellsize = 1.0
        with pytest.raises(ValueError):
            particle_track.unstruct2grid(coords, quantity, cellsize,
                                         k_nearest_neighbors=1)

    def test_no_initialized_particles(self):
        particle = particle_track.Particles(goodparams)
        with pytest.raises(ValueError):
            particle.run_iteration()


def test_coord2ind():
    coords = [(10, 10),
              (10, 15),
              (15, 20)]
    raster_origin = (2, 5)
    raster_size = (25, 25)
    cellsize = 1.0
    inds = particle_track.coord2ind(coords, raster_origin,
                                    raster_size, cellsize)
    assert inds[0] == (20, 8)
    assert inds[1] == (15, 8)
    assert inds[2] == (10, 13)

def test_ind2coord():
    walk_data = dict()
    walk_data['xinds'] = [[20], [15, 10]]
    walk_data['yinds'] = [[8], [8, 13]]
    raster_origin = (2, 5)
    raster_size = (25, 25)
    cellsize = 1.0
    new_data = particle_track.ind2coord(walk_data, raster_origin,
                                        raster_size, cellsize)
    assert new_data['xcoord'][0] == [10.0]
    assert new_data['ycoord'][0] == [10.0]
    assert new_data['xcoord'][1] == [10.0, 15.0]
    assert new_data['ycoord'][1] == [15.0, 20.0]

def test_exposure_time_y():
    walk_data = dict()
    walk_data['xinds'] = [[1, 1, 1, 1, 1]]
    walk_data['yinds'] = [[1, 2, 3, 4, 5]]
    walk_data['travel_times'] = [[2, 4, 6, 8, 10]]
    roi = np.zeros((6, 6))
    roi[:, 2:4] = 1
    exp_times = particle_track.exposure_time(walk_data, roi)
    assert exp_times[0] == 4.0

def test_exposure_time_y():
    walk_data = dict()
    walk_data['yinds'] = [[1, 1, 1, 1, 1]]
    walk_data['xinds'] = [[1, 2, 3, 4, 5]]
    walk_data['travel_times'] = [[2, 4, 6, 8, 10]]
    roi = np.zeros((6, 6))
    roi[2:4, :] = 1
    exp_times = particle_track.exposure_time(walk_data, roi)
    assert exp_times[0] == 4.0

def test_exposure_reenter():
    walk_data = dict()
    walk_data['yinds'] = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
    walk_data['xinds'] = [[1, 2, 3, 4, 5, 4, 3, 2, 3, 4, 5]]
    walk_data['travel_times'] = [[2, 4, 6, 8, 10, 11, 12, 13, 14, 15, 16]]
    roi = np.zeros((6, 6))
    roi[2:3, :] = 1
    exp_times = particle_track.exposure_time(walk_data, roi)
    assert exp_times[0] == 3.0

def test_unstruct2grid_k1():
    coords = [(10.5, 10.1),
              (10.1, 15.1),
              (15.2, 20.2)]
    quantity = [1, 2, 3]
    cellsize = 1.0
    interp_func, gridd = particle_track.unstruct2grid(coords, quantity,
                                                      cellsize,
                                                      k_nearest_neighbors=1)
    assert np.all(gridd == np.array([[3., 3., 3., 3., 3., 3., 3.],
                                     [2., 3., 3., 3., 3., 3., 3.],
                                     [2., 2., 3., 3., 3., 3., 3.],
                                     [2., 2., 2., 3., 3., 3., 3.],
                                     [2., 2., 2., 2., 3., 3., 3.],
                                     [2., 2., 2., 2., 2., 3., 3.],
                                     [2., 2., 2., 2., 2., 2., 3.],
                                     [2., 2., 2., 2., 2., 2., 2.],
                                     [2., 2., 2., 2., 2., 2., 2.],
                                     [1., 1., 1., 1., 1., 1., 1.],
                                     [1., 1., 1., 1., 1., 1., 1.],
                                     [1., 1., 1., 1., 1., 1., 1.]]))

def test_unstruct2grid_k3():
    coords = [(1.5, 1.1),
              (0.1, 2.1),
              (1.2, 2.2)]
    quantity = [1, 2, 3]
    cellsize = 1.0
    interp_func, gridd = particle_track.unstruct2grid(coords, quantity,
                                                      cellsize,
                                                      k_nearest_neighbors=3)
    assert pytest.approx(gridd == np.array([[2.13911589, 2.26676874,
                                             2.1792037],
                                            [2., 2.68254467, 2.10026059],
                                            [1.96968263, 1.6122416,
                                             1.65818041]]))
