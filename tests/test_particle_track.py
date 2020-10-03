from __future__ import division, print_function, absolute_import
from builtins import range, map
from math import floor, sqrt, pi
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
params.seed_xloc = [1]
params.seed_yloc = [1]
params.Np_tracer = 1
params.dx = 1
params.depth = np.ones((3,3))
params.stage = np.ones((3,3))
params.qx = np.zeros((3,3))
params.qy = np.ones((3,3))
params.theta = 1
params.model = 'DeltaRCM'
particle = particle_track.Particles(params)

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
    init_walk_data = particle.generate_particles()
    all_walk_data = particle.run_iteration(init_walk_data)
    assert all_walk_data['xinds'][0][0] == params.seed_xloc[0]

def test_start_pairs_Y():
    init_walk_data = particle.generate_particles()
    all_walk_data = particle.run_iteration(init_walk_data)
    assert all_walk_data['yinds'][0][0] == params.seed_yloc[0]

def test_start_pairs_X2():
    init_walk_data = particle.generate_particles(seed_xloc=params.seed_xloc)
    all_walk_data = particle.run_iteration(init_walk_data)
    assert all_walk_data['xinds'][0][0] == params.seed_xloc[0]

def test_start_pairs_Y2():
    init_walk_data = particle.generate_particles(seed_yloc=params.seed_yloc)
    all_walk_data = particle.run_iteration(init_walk_data)
    assert all_walk_data['yinds'][0][0] == params.seed_yloc[0]

def test_travel_time():
    # Particle doesn't travel in the 3x3 space due to the 'sticky' edge
    # conditions so check that travel time is 0 and particle hasn't moved
    np.random.seed(0)
    init_walk_data = particle.generate_particles()
    all_walk_data = particle.run_iteration(init_walk_data)
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0

def test_init_params():
    # test initialization of params class
    params = particle_track.modelParams()
    # make assertions
    assert params.seed_xloc is None
    assert params.seed_yloc is None
    assert params.Np_tracer is None
    assert params.dx is None
    assert params.depth is None
    assert params.stage is None
    assert params.qx is None
    assert params.qy is None
    assert params.u is None
    assert params.v is None


def test_previous_walk_data():
    # test of loading previously defined walk data
    np.random.seed(0)
    old_init = particle.generate_particles()
    old_walk_data = particle.run_iteration(old_init)
    # try to do another walk - test just makes sure code doesn't break
    np.random.seed(0)
    all_walk_data = particle.run_iteration(init_walk_data=old_walk_data)
    assert all_walk_data['xinds'][0][0] == 1
    assert all_walk_data['yinds'][0][0] == 1
    assert all_walk_data['travel_times'][0][0] == 0.0


class TestValueErrors:
    """
    Catching the ValueErrors and edge cases
    """

    def test_seedx(self):
        params = particle_track.modelParams()
        # expect ValueError
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_seedy(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        # expect ValueError
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_seedx_none(self):
        class pobj():
            pass
        params = pobj
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_seedy_none(self):
        class pobj():
            pass
        params = pobj
        params.seed_xloc = [1]
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_num_tracers(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_dx(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_broken_depth(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.depth = 'badstring'
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_missing_depth(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_depth_via_stage_topo(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
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
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.stage = 'badstring'
        params.topography = np.zeros((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_stage_broken(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.stage = 'badstring'
        params.depth = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_stage_via_depth_topo(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
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
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.topography = 'badstring'
        params.depth = np.ones((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_missing_stage(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.depth = np.ones((3,3))
        params.u = np.ones((3,3))
        params.v = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_rcm_model_uv(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
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
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        params.model = 'DeltaRCM'
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)

    def test_model_q(self):
        params = particle_track.modelParams()
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
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
        params.seed_xloc = [1]
        params.seed_yloc = [1]
        params.Np_tracer = 1
        params.dx = 1
        params.depth = np.ones((3,3))
        params.topography = np.zeros((3,3))
        params.u = np.ones((3,3))
        with pytest.raises(ValueError):
            particle = particle_track.Particles(params)


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
