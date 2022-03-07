from __future__ import division, print_function, absolute_import
from builtins import range, map
from math import floor, sqrt, pi
import pytest
import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np
import dorado.particle_track as pt
from dorado.particle_track import Particles


class TestRCM:
    """Class to test particles on RCM field."""

    # set up DeltaRCM test case
    rcm_data = np.load('tests/data/ex_deltarcm_data.npz')
    params = pt.modelParams()
    params.stage = rcm_data['stage']
    params.depth = rcm_data['depth']
    params.dx = 50.
    # just do one particle in a set location for this test
    seed_xloc = [16]
    seed_yloc = [140]
    Np_tracer = 1
    # no discharge data so use 0s as substitute
    params.qx = np.zeros_like(params.depth)
    params.qy = np.zeros_like(params.depth)
    params.model = 'DeltaRCM'

    def test_few_steps_RCM(self):
        '''
        Test running a few steps
        '''
        # defining / initializing
        particle = Particles(self.params)  # define the particle
        np.random.seed(0)  # fix the random seed
        # init the walk data
        particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)

        # 3 iterations
        for i in list(range(0, 3)):
            all_walk_data = particle.run_iteration()

        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # check all positions and times
        assert all_walk_data['xinds'][0] == [16, 17, 18, 19]
        assert all_walk_data['yinds'][0] == [140, 141, 141, 141]
        assert all_walk_data['travel_times'][0] == \
            [0, 7007593448.337235, 10439964733.462337, 13698473384.876724]

    def test_set_time_RCM_previousdata(self):
        '''
        Test setting a time target when using old walk data
        '''
        # defining / initializing
        particle = Particles(self.params)  # define the particle
        np.random.seed(0)  # fix the random seed

        # generate the particles
        particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)  # init walk data

        # set time
        all_walk_data = particle.run_iteration(target_time=5e6)
        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # check all positions and times
        assert all_walk_data['xinds'][0] == [16, 17]
        assert all_walk_data['yinds'][0] == [140, 141]
        assert all_walk_data['travel_times'][0] == [0, 7007593448.337235]

    def test_set_time_RCM(self):
        '''
        Test setting a time target
        '''
        # defining / initializing
        particle = Particles(self.params)  # define the particle
        np.random.seed(0)  # fix the random seed

        # generate particles
        particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)  # init walk data

        # set time
        all_walk_data = particle.run_iteration(target_time=5e6)

        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # check all positions and times
        assert all_walk_data['xinds'][0] == [16, 17]
        assert all_walk_data['yinds'][0] == [140, 141]
        assert all_walk_data['travel_times'][0] == [0, 7007593448.337235]


class TestANUGA:
    """Class to test particles on ANUGA flow field."""

    # set up anuga test case - going to use subset of domain
    an_data = np.load('tests/data/ex_anuga_data.npz')
    an_params = pt.modelParams()
    an_params.stage = an_data['depth'][40:61, 40:61]  # just using depth here
    an_params.depth = an_data['depth'][40:61, 40:61]
    # just do one particle in a set location for this test
    seed_xloc = [5]
    seed_yloc = [10]
    Np_tracer = 1
    an_params.dx = 50.
    # no discharge data so use 0s as substitute
    an_params.qx = an_data['qx'][40:61, 40:61]
    an_params.qy = an_data['qy'][40:61, 40:61]
    an_params.model = 'Anuga'

    def test_few_steps_anuga(self):
        '''
        Test running a few steps
        '''
        # defining / initializing
        an_particle = Particles(self.an_params)  # define the particle
        np.random.seed(0)  # fix the random seed
        an_particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)  # init walk data

        # 3 iterations
        for i in list(range(0, 3)):
            all_walk_data = an_particle.run_iteration()

        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # check all positions and times
        assert all_walk_data['xinds'][0] == [5, 6, 7, 8]
        assert all_walk_data['yinds'][0] == [10, 10, 10, 10]
        assert all_walk_data['travel_times'][0] == \
            [0, 132.8105258698104, 258.1975289086354, 375.18128832888203]

    def test_boundary_anuga(self):
        '''
        Test running into the boundary
        '''
        # defining / initializing
        an_particle = Particles(self.an_params)  # define the particle
        np.random.seed(0)  # fix the random seed
        an_particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)  # init walk data

        # 20 iterations
        for i in list(range(0, 20)):
            all_walk_data = an_particle.run_iteration()

        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # particle reaches boundary after 18 iterations
        # data shouldn't be recorded after boundary is reached
        assert len(all_walk_data['xinds'][0]) == 15
        assert len(all_walk_data['yinds'][0]) == 15
        assert len(all_walk_data['travel_times'][0]) == 15

    def test_boundary_travel_time_anuga(self):
        '''
        Test running into the boundary and not reaching travel time target
        '''
        # defining / initializing
        an_particle = Particles(self.an_params)  # define the particle
        np.random.seed(0)  # fix the random seed

        # generate particles
        an_particle.generate_particles(
            self.Np_tracer, self.seed_xloc, self.seed_yloc)  # init walk data

        # set target time for iterations
        all_walk_data = an_particle.run_iteration(target_time=1000.0)

        # make assertions
        # check initial position and travel time
        assert all_walk_data['xinds'][0][0] == self.seed_xloc[0]
        assert all_walk_data['yinds'][0][0] == self.seed_yloc[0]
        assert all_walk_data['travel_times'][0][0] == 0.0
        # particle reaches boundary after 18 iterations
        # data shouldn't be recorded after boundary is reached
        assert len(all_walk_data['xinds'][0]) == 9
        assert len(all_walk_data['yinds'][0]) == 9
        assert len(all_walk_data['travel_times'][0]) == 9
