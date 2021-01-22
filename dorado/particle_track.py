# -*- coding: utf-8 -*-
"""
Lower-level methods to manage parameters and particle movement.

Particle class for managing the definition of particle attributes and
parameters of the domain as well as iterative movement of the particles
through the domain.

Project Homepage: https://github.com/passaH2O/dorado
"""
from __future__ import division, print_function, absolute_import
from builtins import range
from math import pi
import numpy as np
from tqdm import tqdm
import dorado.lagrangian_walker as lw


class modelParams:
    """Parameter class with attributes and grid particles will be routed on.

    The parameters class, :obj:`dorado.particle_track.modelParams`,
    must be populated with user-defined
    attributes of the grid the particles will be modeled on.

    **Required Parameters:**

        dx : `float`
            Length along one square cell face

        depth : `numpy.ndarray`
            Array of water depth values, if absent then the stage and
            topography arrays will be used to compute it

        stage : `numpy.ndarray`
            Array of water stage values, if absent then the depth and
            topography arrays will be used to compute it

        qx : `numpy.ndarray`
            Array of the x-component of flow discharge

        qy : `numpy.ndarray`
            Array of the y-component of flow discharge

        u : `numpy.ndarray`
            Array of the x-component of flow velocity

        v : `numpy.ndarray`
            Array of the y-component of flow velocity

    **Optional Parameters:**

        topography : `numpy.ndarray`
            Array of cell elevation values

        model : `str`
            Name of the hydrodynamic model input being used (e.g. 'DeltaRCM')

        theta : `float`
            First of two weighting parameters for the weighted random walk.
            Default value is 1.0, higher values give higher weighting
            probabilities to cells with greater water depths

        gamma  : `float`
            Second of two weighting parameters for the weighted random walk.
            Default value is 0.05. Gamma must be in the range [0,1]. Gamma == 0
            means that the random walk weights are independent of the discharge
            values, and instead are based on the water surface gradient (the
            stage). Gamma == 1 means that the random walk weights are not
            dependent on the surface gradient, and instead are based on the
            inertial forces (the flow discharge).

        diff_coeff : `float`
            Diffusion/dispersion coefficient for use in travel time
            computation. If set to 0.0, flow is purely advection with no
            diffusion. Higher values lead to more spread in exit age
            distribution. Max diffusion time in any given step is
            0.5*diff_coeff percent. Default value is 0.2 (i.e. max of 10%)

        dry_depth : `float`
            Minimum depth for a cell to be considered wet,
            default value is 0.1m

        cell_type : `numpy.ndarray`
            Array of the different types of cells in the domain where 2 = land,
            1 = channel, 0 = ocean, and -1 = edge. If not explicitly defined
            then the values are estimated based on the depth array and the
            defined dry_depth

        steepest_descent : `bool`
            Toggle for routing based on a steepest descent rather than the
            weighted random walk. If True, then the highest weighted cells are
            used to route the particles. Default value is False.

    This list of expected parameter values can also be obtained by querying the
    class attributes with `dir(modelParams)`, `modelParams.__dict__`, or
    `vars(modelParams)`.

    """

    def __init__(self):
        """Create the expected variables for the modelParams class.

        Variables are initialized as NoneType they need to be assigned by the
        user. Due to the wide variety of data formats produced by differeny
        hydrodynamic models, this is process is not automated and must be
        handled on a case-by-case basis.

        """
        self.dx = None
        self.depth = None
        self.stage = None
        self.qx = None
        self.qy = None
        self.u = None
        self.v = None


class Particles():
    """Class for the particle(s) that is(are) going to be routed."""

    def __init__(self, params):
        """Check input parameters and assign default values where/if needed.

        Methods require a class of parameters
        (:obj:`dorado.particle_track.modelParams`) to be passed
        to the Particles class. e.g. particle = Particles(modelParams)

        Initialization tries to assign each value from the parameter class,
        otherwise an error is raised or default values are assigned when
        possible/sensible

        """
        # REQUIRED PARAMETERS #
        # Define the length along one cell face (assuming square cells)
        if getattr(params, 'dx', None) is None:
            raise ValueError("Length of cell face (modelParams.dx) is"
                             " undefined")
        else:
            self.dx = params.dx

        # Define the water depth array
        if getattr(params, 'depth', None) is not None:
            try:
                self.depth = params.depth
                self.depth[np.isnan(self.depth)] = 0
            except Exception:
                raise ValueError("Water depth array incorrectly defined.")
        elif getattr(params,
                     'stage',
                     None) is not None and getattr(params,
                                                   'topography',
                                                   None) is not None:
            try:
                self.depth = params.stage - params.topography
                self.depth[self.depth < 0] = 0
                self.depth[np.isnan(self.depth)] = 0
            except Exception:
                raise ValueError("Insufficient information: Specify depth")
        else:
            raise ValueError("Insufficient information: Specify depth")

        # Define the water stage array
        if getattr(params, 'stage', None) is not None:
            try:
                self.stage = params.stage
                self.stage[self.depth == 0] = np.nan
            except Exception:
                raise ValueError("Water stage array incorrectly defined.")
        elif getattr(params,
                     'topography',
                     None) is not None and getattr(params,
                                                   'depth',
                                                   None) is not None:
            try:
                self.stage = params.topography + params.depth
                self.stage[self.depth == 0] = np.nan
            except Exception:
                raise ValueError("Insufficient information: Specify stage")
        else:
            raise ValueError("Insufficient information: Specify stage")

        # check if hydrodynamic model input has been specified
        if getattr(params, 'model', None) is not None:
            pass
        else:
            params.model = []

        # Define discharge and velocities for all cells in domain
        if params.model == 'DeltaRCM':
            if params.qx is not None and params.qy is not None:
                self.qx = params.qx
                self.qx[np.isnan(self.qx)] = 0
                self.u = self.qx*self.depth/(self.depth**2 + 1e-8)
                self.u[np.isnan(self.u)] = 0

                self.qy = params.qy
                self.qy[np.isnan(self.qy)] = 0
                self.v = self.qy*self.depth/(self.depth**2 + 1e-8)
                self.v[np.isnan(self.v)] = 0

            elif params.u is not None and params.v is not None:
                self.u = params.u
                self.u[np.isnan(self.u)] = 0
                self.qx = self.u*self.depth

                self.v = params.v
                self.v[np.isnan(self.v)] = 0
                self.qy = self.v*self.depth

            else:
                raise ValueError("Insufficient information:"
                                 " Specify velocities/discharge")
        else:
            if params.qx is not None and params.qy is not None:
                self.qx = -1*params.qy
                self.qx[np.isnan(self.qx)] = 0
                self.u = self.qx*self.depth/(self.depth**2 + 1e-8)
                self.u[np.isnan(self.u)] = 0

                self.qy = params.qx
                self.qy[np.isnan(self.qy)] = 0
                self.v = self.qy*self.depth/(self.depth**2 + 1e-8)
                self.v[np.isnan(self.v)] = 0

            elif params.u is not None and params.v is not None:
                self.u = -1*params.v
                self.u[np.isnan(self.u)] = 0
                self.qx = self.u*self.depth

                self.v = params.u
                self.v[np.isnan(self.v)] = 0
                self.qy = self.v*self.depth

            else:
                raise ValueError("Insufficient information:"
                                 " Specify velocities/discharge")

        # Define field of velocity magnitude (for travel time calculation)
        self.velocity = np.sqrt(self.u**2+self.v**2)
        # cannot have 0/nans - leads to infinite/nantravel times
        self.velocity[self.velocity < 1e-8] = 1e-8
        self.u[np.abs(self.u) < 1e-8] = 1e-8
        self.v[np.abs(self.v) < 1e-8] = 1e-8
        # Compute velocity orientation at each cell
        self.velocity_angle = np.arctan2(-1.0*self.u, self.v)

        # OPTIONAL PARAMETERS (Have default values) #
        # Define the theta used to weight the random walk
        # Higher values give higher weighting probabilities to deeper cells
        try:
            self.theta = float(params.theta)
        except Exception:
            print("Theta parameter not specified - using 1.0")
            self.theta = 1.0  # if unspecified use 1

        # Gamma parameter used to weight the random walk
        # Sets weight ratio (between 0 and 1):
        # 1 = water surface gradient only (stage based)
        # 0 = inertial force only (discharge based)
        try:
            self.gamma = float(params.gamma)
        except Exception:
            print("Gamma parameter not specified - using 0.05")
            self.gamma = 0.05

        try:
            if params.diff_coeff < 0:
                print("Warning: Specified diffusion coefficient is negative."
                      " Rounding up to zero")
                params.diff_coeff = 0.0
            elif params.diff_coeff >= 2:
                print("Warning: Diffusion behaves non-physically when"
                      " coefficient >= 2")
            self.diff_coeff = float(params.diff_coeff)
        except Exception:
            if getattr(params, 'steepest_descent', False) is True:
                print("Diffusion disabled for steepest descent")
                self.diff_coeff = 0.0
            else:
                print("Diffusion coefficient not specified - using 0.2")
                self.diff_coeff = 0.2

        # Minimum depth for cell to be considered wet
        try:
            self.dry_depth = params.dry_depth
        except Exception:
            print("minimum depth for wetness not defined - using 10 cm")
            self.dry_depth = 0.1

        # Cell types: 2 = land, 1 = channel, 0 = ocean, -1 = edge
        try:
            self.cell_type = params.cell_type
        except Exception:
            print("Cell Types not specified - Estimating from depth")
            self.cell_type = np.zeros_like(self.depth, dtype='int')
            self.cell_type[self.depth < self.dry_depth] = 2
            self.cell_type = np.pad(self.cell_type[1:-1, 1:-1], 1, 'constant',
                                    constant_values=-1)

        # Steepest descent toggle - turns off randomness and uses highest
        # weighted value instead of doing weighted random walk
        # note: chooses randomly in event of ties
        try:
            if params.steepest_descent is True:
                print("Using steepest descent")
                self.steepest_descent = True
            else:
                print("Using weighted random walk")
                self.steepest_descent = False
        except Exception:
            print("Using weighted random walk")
            self.steepest_descent = False

        # DEFAULT PARAMETERS (Can be defined otherwise) #

        sqrt2 = np.sqrt(2)
        sqrt05 = np.sqrt(0.5)

        # Define distances between cells in D8 sense
        try:
            self.distances = params.distances
        except Exception:
            # defined if not given
            self.distances = np.array([[sqrt2, 1, sqrt2],
                                       [1, 1, 1],
                                       [sqrt2, 1, sqrt2]])

        # D8 components of x-unit vector
        try:
            self.ivec = params.ivec
        except Exception:
            # defined if not given
            self.ivec = np.array([[-sqrt05, 0, sqrt05],
                                  [-1, 0, 1],
                                  [-sqrt05, 0, sqrt05]])

        # D8 components of y-unit vector
        try:
            self.jvec = params.jvec
        except Exception:
            # defined if not given
            self.jvec = np.array([[-sqrt05, -1, -sqrt05],
                                  [0, 0, 0],
                                  [sqrt05, 1, sqrt05]])

        # Positive/Negative x-directions
        try:
            self.iwalk = params.iwalk
        except Exception:
            # defined if not given
            self.iwalk = np.array([[-1, 0, 1],
                                   [-1, 0, 1],
                                   [-1, 0, 1]])

        # Positive/Negative y-directions
        try:
            self.jwalk = params.jwalk
        except Exception:
            # defined if not given
            self.jwalk = np.array([[-1, -1, -1],
                                   [0, 0, 0],
                                   [1, 1, 1]])

        # Angles of D8 step directions
        try:
            self.angles = params.angles
        except Exception:
            # defined if not given
            self.angles = np.array([[3*pi/4, pi/2, pi/4],
                                    [pi, 0, 0],
                                    [5*pi/4, 3*pi/2, 7*pi/4]])

        # initialize number of particles as 0
        self.Np_tracer = 0

        # initialize the walk_data
        self.walk_data = None

        # initialize routing weights array
        self.weight = np.zeros((self.stage.shape[0], self.stage.shape[1], 9))

    # function to clear walk data if you've made a mistake while generating it
    def clear_walk_data(self):
        """Manually reset self.walk_data back to None."""
        self.walk_data = None

    # generate a group of particles in a set of initial locations
    def generate_particles(self,
                           Np_tracer,
                           seed_xloc,
                           seed_yloc,
                           method='random',
                           previous_walk_data=None):
        """Generate a set of particles in defined locations.

        After a :obj:`dorado.particle_track.Particles` class has been
        initialized, this function can be called to create some particles
        within a given region. If particles are to be seeded in different
        groupings in different regions, this function can be called multiple
        times in succession prior to moving them via
        :obj:`dorado.particle_track.run_iteration`.

        Walk data is stored within :obj:`dorado.particle_track.Particles` as
        `Particles.walk_data` so if this method is called in succession, even
        without the `previous_walk_data` flag, the new particle seed locations
        will be appended to the `walk_data` attribute of
        :obj:`dorado.particle_track.Particles`.

        **Inputs** :

            Np_tracer : `int`
                Number of particles to generate.

            seed_xloc : `list`
                List of x-coordinates over which to initially distribute the
                particles.

            seed_yloc : `list`
                List of y-coordinates over which to initially distribute the
                particles.

            method : `str`, optional
                Specify the type of particle generation you wish to use.
                Current options are 'random' and 'exact' for seeding particles
                either randomly across a set of points, or exactly at a set
                of defined locations. If unspecified, the default is 'random'.

            previous_walk_data : `dict`, optional
                Dictionary of all prior x locations, y locations, and travel
                times. This input parameter should only be used if 'new' or
                'indpendent' walk data exists (e.g. data loaded from a .json
                file). Otherwise the walk data stored in `self.walk_data` is
                likely to contain the previous information. Order of indices is
                previous_walk_data[field][particle][iter], where e.g.
                ['travel_times'][5][10] is the travel time of the 5th particle
                at the 10th iteration.

        """
        # if the values in self are invalid the input checked will catch them
        # do input type checking
        Np_tracer, seed_xloc, seed_yloc = gen_input_check(Np_tracer,
                                                          seed_xloc,
                                                          seed_yloc,
                                                          method)

        init_walk_data = dict()  # create init_walk_data dictionary

        # initialize new travel times list
        new_start_times = [0.]*Np_tracer

        if method == 'random':
            # create random start indices based on x, y locations and np_tracer
            new_start_xindices = [lw.random_pick_seed(seed_xloc) for x
                                  in list(range(Np_tracer))]
            new_start_yindices = [lw.random_pick_seed(seed_yloc) for x
                                  in list(range(Np_tracer))]

        elif method == 'exact':
            # create exact start indices based on x, y locations and np_tracer
            # get number of particles per location and remainder left out
            Np_divy = divmod(Np_tracer, len(seed_xloc))
            # assign the start locations for the evenly divided particles
            new_start_xindices = seed_xloc * Np_divy[0]
            new_start_yindices = seed_yloc * Np_divy[0]
            # add the remainder ones to the list one by one via looping
            for i in range(0, Np_divy[1]):
                new_start_xindices = new_start_xindices + [seed_xloc[i]]
                new_start_yindices = new_start_yindices + [seed_yloc[i]]

        # Now initialize vectors that will create the structured list
        new_xinds = [[new_start_xindices[i]] for i in
                     list(range(Np_tracer))]
        new_yinds = [[new_start_yindices[i]] for i in
                     list(range(Np_tracer))]
        new_times = [[new_start_times[i]] for i in list(range(Np_tracer))]

        # establish the start values
        start_xindices = new_xinds
        start_yindices = new_yinds
        start_times = new_times

        if self.walk_data is not None:
            # if there is walk_data from a previous call to the generator,
            # or from simulating particle transport previously, then we want
            # to keep it and append the new data to it
            internal_xinds = self.walk_data['xinds']
            internal_yinds = self.walk_data['yinds']
            internal_times = self.walk_data['travel_times']

            # combine internal and new lists of particle information
            start_xindices = internal_xinds + start_xindices
            start_yindices = internal_yinds + start_yindices
            start_times = internal_times + start_times

        if previous_walk_data is not None:
            # If the generator has been run before, or if new
            # particles are going to be added to a set of pre-existing
            # particles, then the walk_data dictionary from either can
            # be fed into this function as input
            prev_xinds = previous_walk_data['xinds']
            prev_yinds = previous_walk_data['yinds']
            prev_times = previous_walk_data['travel_times']

            # combine previous and new lists of particle information
            start_xindices = prev_xinds + start_xindices
            start_yindices = prev_yinds + start_yindices
            start_times = prev_times + start_times

        # determine the new total number of particles we have now
        self.Np_tracer = len(start_xindices)

        # store information in the init_walk_data dictionary and return it
        init_walk_data['xinds'] = start_xindices
        init_walk_data['yinds'] = start_yindices
        init_walk_data['travel_times'] = start_times

        # store the initialized walk data within self
        self.walk_data = init_walk_data

    # run an iteration where particles are moved
    # have option of specifying the particle start locations
    # otherwise they are randomly placed within x and y seed locations
    def run_iteration(self,
                      target_time=None):
        """Run an iteration of the particle routing.

        Runs an iteration of the particle routing.
        Returns at each step the particle's locations and travel times.

        **Inputs** :

            target_time : `float`, optional
                The travel time (seconds) each particle should aim to have at
                end of this iteration. If left undefined, then just one
                iteration is run and the particles will be out of sync in time.
                Note that this loop will terminate before the target_time if
                the particle exceeds the hard-coded limit of 1e4 steps

        **Outputs** :

            all_walk_data : `dict`
                Dictionary of all x and y locations and travel times. Order of
                indices is previous_walk_data[field][particle][iter],
                where ['travel_times'][5][10] is the travel time of the 5th
                particle at the 10th iteration

        """
        # if there is no walk data raise an error
        if self.walk_data is None:
            raise ValueError('Particles have not been initialized.'
                             ' Call the `particle_generator()` function.')

        all_walk_data = dict()  # init all_walk_data dictionary

        # if self.Np_tracer is not already defined, we can define it from
        # the init_walk_data
        if self.Np_tracer == 0:
            self.Np_tracer = len(self.walk_data['xinds'])

        # read in information from the init_walk_data dictionary
        # most recent locations
        all_xinds = self.walk_data['xinds']
        start_xindices = [all_xinds[i][-1] for i in
                          list(range(self.Np_tracer))]
        all_yinds = self.walk_data['yinds']
        start_yindices = [all_yinds[i][-1] for i in
                          list(range(self.Np_tracer))]
        # most recent times
        all_times = self.walk_data['travel_times']
        start_times = [all_times[i][-1] for i in
                       list(range(self.Np_tracer))]

        # merge x and y indices into list of [x,y] pairs
        start_pairs = [[start_xindices[i], start_yindices[i]] for i in
                       list(range(self.Np_tracer))]

        # Do the particle movement
        if target_time is None:
            # If we're not aiming for a specific time, step the particles
            new_inds, travel_times = lw.particle_stepper(self, start_pairs,
                                                         start_times)

            for ii in list(range(self.Np_tracer)):
                # Don't duplicate location
                # if particle is standing still at a boundary
                if new_inds[ii] != start_pairs[ii]:
                    # Append new information
                    all_xinds[ii].append(new_inds[ii][0])
                    all_yinds[ii].append(new_inds[ii][1])
                    all_times[ii].append(travel_times[ii])

            # Store travel information in all_walk_data
            all_walk_data['xinds'] = all_xinds
            all_walk_data['yinds'] = all_yinds
            all_walk_data['travel_times'] = all_times

        else:
            # If we ARE aiming for a specific time
            # iterate each particle until we get there
            # Loop through all particles
            for ii in list(range(self.Np_tracer)):
                if len(all_times[ii]) > 1:
                    est_next_dt = all_times[ii][-1] - all_times[ii][-2]
                else:
                    # Initialize a guess for the next iteration's timestep
                    est_next_dt = 0.1
                count = 1

                # Only iterate if this particle isn't already at a boundary:
                if -1 not in self.cell_type[all_xinds[ii][-1]-1:
                                            all_xinds[ii][-1]+2,
                                            all_yinds[ii][-1]-1:
                                            all_yinds[ii][-1]+2]:
                    # Loop until |target time - current time| <
                    #            |target time - estimated next time|
                    while abs(all_times[ii][-1] - target_time) >= \
                          abs(all_times[ii][-1] + est_next_dt - target_time):
                        # for particle ii, take a step from most recent index
                        new_inds, travel_times = lw.particle_stepper(self,
                                                    [[all_xinds[ii][-1],
                                                      all_yinds[ii][-1]]],
                                                    [all_times[ii][-1]])

                        # Don't duplicate location
                        # if particle is standing still at a boundary
                        if new_inds[0] != [all_xinds[ii][-1],
                                           all_yinds[ii][-1]]:
                            all_xinds[ii].append(new_inds[0][0])
                            all_yinds[ii].append(new_inds[0][1])
                            all_times[ii].append(travel_times[0])
                        else:
                            break

                        # Use that timestep to estimate how long the
                        # next one will take
                        est_next_dt = max(0.1, all_times[ii][-1] -
                                          all_times[ii][-2])
                        count += 1
                        if count > 1e4:
                            print('Warning: Particle iterations exceeded limit'
                                  ' before reaching target time. Try smaller'
                                  ' time-step')
                            break

            # Store travel information in all_walk_data
            all_walk_data['xinds'] = all_xinds
            all_walk_data['yinds'] = all_yinds
            all_walk_data['travel_times'] = all_times

            # re-write the walk_data attribute of self
            self.walk_data = all_walk_data

        return all_walk_data


def gen_input_check(Np_tracer, seed_xloc, seed_yloc, method):
    """Check the inputs provided to :obj:`generate_particles()`.

    This function does input type checking and either succeeds or returns
    an error if the types provided cannot be converted to those that we
    expect in :obj:`dorado.particle_track.generate_particles()`.

    **Inputs** :

        Np_tracer : `int`
            Number of particles to generate.

        seed_xloc : `list`
            List of x-coordinates over which to initially distribute the
            particles

        seed_yloc : `list`
            List of y-coordinates over which to initially distribute the
            particles

        method : `str`, optional
            Type of particle generation to use, either 'random' or 'exact'.

    **Outputs** :

        Np_tracer : `int`
            Number of particles to generate.

        seed_xloc : `list`
            List of x-coordinates over which to initially distribute the
            particles

        seed_yloc : `list`
            List of y-coordinates over which to initially distribute the
            particles

    """
    try:
        Np_tracer = int(Np_tracer)
    except Exception:
        raise TypeError("Np_tracer input type was not int.")
    try:
        # try to flatten if this is an array that is multi-dimensional
        if len(np.shape(seed_xloc)) > 1:
            seed_xloc = seed_xloc.flatten()
        # try to coerce into a list
        seed_xloc = [int(x) for x in list(seed_xloc)]
    except Exception:
        raise TypeError("seed_xloc input type was not a list.")
    try:
        # try to flatten if this is an array that is multi-dimensional
        if len(np.shape(seed_yloc)) > 1:
            seed_yloc = seed_yloc.flatten()
        # try to coerce into a list
        seed_yloc = [int(y) for y in list(seed_yloc)]
    except Exception:
        raise TypeError("seed_yloc input type was not a list.")
    if isinstance(method, str) is False:
        raise TypeError('Method provided was not a string.')
    elif method not in ['random', 'exact']:
        raise ValueError('Method input is not a valid method, must be'
                         ' "random" or "exact".')

    return Np_tracer, seed_xloc, seed_yloc


def coord2ind(coordinates, raster_origin, raster_size, cellsize):
    """Convert geographical coordinates into raster index coordinates.

    Assumes geographic coordinates are projected onto a Cartesian grid.
    Accepts coordinates in meters or decimal degrees.

    **Inputs** :

        coordinates : `list`
            List [] of (x,y) pairs or tuples of coordinates to be
            converted from starting units (e.g. meters UTM) into raster
            index coordinates used in particle routing functions.

        raster_origin : `tuple`
            Tuple of the (x,y) raster origin in physical space, i.e. the
            coordinates of lower left corner. For rasters loaded from a
            GeoTIFF, lower left corner can be obtained using e.g. gdalinfo

        raster_size : `tuple`
            Tuple (L,W) of the raster dimensions, i.e. the output of
            numpy.shape(raster).

        cellsize : `float or int`
            Length along one square cell face.

    **Outputs** :

        inds : `list`
            List [] of tuples (x,y) of raster index coordinates

    """
    x_orig = float(raster_origin[0])
    y_orig = float(raster_origin[1])

    cellsize = float(cellsize)

    L = int(raster_size[0])  # Need domain extent

    inds = []
    for i in list(range(0, len(coordinates))):
        # Do coordinate transform:
        new_ind = (int(L - round((coordinates[i][1] - y_orig)/cellsize)),
                   int(round((coordinates[i][0] - x_orig)/cellsize)))
        inds.append(new_ind)

    return inds


def ind2coord(walk_data, raster_origin, raster_size, cellsize):
    """Convert raster index coordinates into geographical coordinates.

    Appends the walk_data dictionary from the output of run_iteration
    with additional fields 'xcoord' and 'ycoord' in projected geographic
    coordinate space. Locations align with cell centroids.
    Units of output coordinates match those of raster_origin and cellsize,
    can be meters or decimal degrees.

    **Inputs** :

        walk_data : `dict`
            Dictionary of all prior x locations, y locations, and travel
            times (the output of run_iteration)

        raster_origin : `tuple`
            Tuple of the (x,y) raster origin in physical space, i.e. the
            coordinates of lower left corner. For rasters loaded from a
            GeoTIFF, lower left corner can be obtained using e.g. gdalinfo

        raster_size : `tuple`
            Tuple (L,W) of the raster dimensions, i.e. the output of
            numpy.shape(raster).

        cellsize : `float or int`
            Length along one square cell face.

    **Outputs** :

        walk_data : `dict`
            Same as the input walk_data dictionary, with the added
            'xcoord' and 'ycoord' fields representing the particles
            geographic position at each iteration.

    """
    x_orig = float(raster_origin[0])
    y_orig = float(raster_origin[1])

    cellsize = float(cellsize)

    L = int(raster_size[0])  # Need domain extent
    Np_tracer = len(walk_data['xinds'])  # Get number of particles

    all_xcoord = []  # Initialize
    all_ycoord = []

    for i in list(range(0, Np_tracer)):
        # Do coordinate transform:
        this_ycoord = [(L-float(j))*cellsize+y_orig for j in
                       walk_data['xinds'][i]]
        all_ycoord.append(this_ycoord)

        this_xcoord = [float(j)*cellsize+x_orig for j in walk_data['yinds'][i]]
        all_xcoord.append(this_xcoord)

    # Save back into dict:
    walk_data['xcoord'] = all_xcoord
    walk_data['ycoord'] = all_ycoord

    return walk_data


def exposure_time(walk_data,
                  region_of_interest):
    """Measure exposure time distribution of particles in a specified region.

    Function to measure the exposure time distribution (ETD) of particles to
    the specified region. For steady flows, the ETD is exactly equivalent to
    the residence time distribution. For unsteady flows, if particles make
    multiple excursions into the region, all of those times are counted.

    **Inputs** :

        walk_data : `dict`
            Output of a previous function call to run_iteration.

        region_of_interest : `int array`
            Binary array the same size as input arrays in modelParams class
            with 1's everywhere inside the region in which we want to
            measure exposure time, and 0's everywhere else.

    **Outputs** :

        exposure_times : `list`
            List of exposure times to region of interest, listed
            in order of particle ID.

    """
    # Initialize arrays to record exposure time of each particle
    Np_tracer = len(walk_data['xinds'])  # Number of particles
    # Array to be populated
    exposure_times = np.zeros([Np_tracer], dtype='float')

    # Loop through particles to measure exposure time
    for ii in tqdm(list(range(0, Np_tracer)), ascii=True):
        # Determine the starting region for particle ii
        previous_reg = region_of_interest[int(walk_data['xinds'][ii][0]),
                                          int(walk_data['yinds'][ii][0])]

        # Loop through iterations
        for jj in list(range(1, len(walk_data['travel_times'][ii]))):
            # Determine the new region and compare to previous region
            current_reg = region_of_interest[int(walk_data['xinds'][ii][jj]),
                                             int(walk_data['yinds'][ii][jj])]

            # Check to see if whole step was inside ROI
            # If so, travel time of the whole step added to ET
            if (current_reg + previous_reg) == 2:
                exposure_times[ii] += (walk_data['travel_times'][ii][jj] -
                                       walk_data['travel_times'][ii][jj-1])
            # Check to see if half of the step was inside ROI
            # (either entering or exiting)
            # If so, travel time of half of the step added to ET
            elif (current_reg + previous_reg) == 1:
                exposure_times[ii] += 0.5*(walk_data['travel_times'][ii][jj] -
                                           walk_data['travel_times'][ii][jj-1])

            # Update previous region
            previous_reg = current_reg

            # Check if particle is still stuck in ROI at the end of the run
            # (which can bias result)
            if jj == len(walk_data['travel_times'][ii])-1:
                if current_reg == 1:
                    print('Warning: Particle ' + str(ii) + ' is still within'
                          ' ROI at final timestep. \n' +
                          'Run more iterations to get tail of ETD')

    return exposure_times.tolist()


def nourishment_area(walk_data, raster_size, sigma=0.7, clip=99.5):
    """Determine the nourishment area of a particle injection

    Function will measure the regions of the domain 'fed' by a seed location,
    as indicated by the history of particle travel locations in walk_data.
    Returns a heatmap raster, in which values indicate number of occasions
    each cell was occupied by a particle (after spatial filtering to reduce
    noise in the random walk, and normalizing by number of particles).

    **Inputs** :

        walk_data : `dict`
            Dictionary of all prior x locations, y locations, and travel
            times (the output of run_iteration)

        raster_size : `tuple`
            Tuple (L,W) of the domain dimensions, i.e. the output of
            numpy.shape(raster).

        sigma : `float`, optional
            Degree of spatial smoothing of the area, implemented using a
            Gaussian kernal of the same sigma, via scipy.ndimage.gaussian_filter
            Default is light smoothing with sigma = 0.7 (to turn off smoothing,
            set sigma = 0)

        clip : `float`, optional
            Percentile at which to truncate the distribution. Particles which
            get stuck can lead to errors at the high-extreme, so this 
            parameter is used to normalize by a "near-max". Default is the
            99.5th percentile. To use true max, specify clip = 100

    **Outputs** :

        visit_freq : `numpy.ndarray`
            Array of normalized particle visit frequency, with cells in the
            range [0, 1] representing the number of instances particles visited
            that cell. If sigma > 0, the array values include spatial filtering

    """
    if sigma > 0:
        from scipy.ndimage import gaussian_filter
    
    # Measure visit frequency
    visit_freq = np.zeros(raster_size)
    for ii in list(range(len(walk_data['xinds']))):
        for jj in list(range(len(walk_data['xinds'][ii]))):
            visit_freq[walk_data['xinds'][ii][jj],
                       walk_data['yinds'][ii][jj]] += 1

    # Clip out highest percentile to correct possible outliers
    vmax = float(np.nanpercentile(visit_freq, clip))
    visit_freq = np.clip(visit_freq, 0, vmax)

    # If applicable, do smoothing
    if sigma > 0:
        visit_freq = gaussian_filter(visit_freq, sigma=sigma)
        visit_freq[visit_freq==np.min(visit_freq)] = np.nan
    else:
        visit_freq[visit_freq==0] = np.nan

    # Normalize to 0-1
    visit_freq = visit_freq/np.nanmax(visit_freq)

    return visit_freq


def nourishment_time(walk_data, raster_size, sigma=0.7, clip=99.5):
    """Determine the nourishment time of a particle injection

    Function will measure the average length of time particles spend in each
    area of the domain for a given seed location, as indicated by the history
    of particle travel times in walk_data. Returns a heatmap raster, in
    which values indicate the average length of time particles tend to remain
    in each cell (after spatial filtering to reduce noise in the random walk).

    **Inputs** :

        walk_data : `dict`
            Dictionary of all prior x locations, y locations, and travel
            times (the output of run_iteration)

        raster_size : `tuple`
            Tuple (L,W) of the domain dimensions, i.e. the output of
            numpy.shape(raster).

        sigma : `float`, optional
            Degree of spatial smoothing of the area, implemented using a
            Gaussian kernal of the same sigma, via scipy.ndimage.gaussian_filter
            Default is light smoothing with sigma = 0.7 (to turn off smoothing,
            set sigma = 0)

        clip : `float`, optional
            Percentile at which to truncate the distribution. Particles which
            get stuck can lead to errors at the high-extreme, so this 
            parameter is used to normalize by a "near-max". Default is the
            99.5th percentile. To use true max, specify clip = 100

    **Outputs** :

        mean_time : `numpy.ndarray`
            Array of mean occupation time, with cell values representing the 
            mean time particles spent in that cell. If sigma > 0, the array
            values include spatial filtering

    """
    if sigma > 0:
        from scipy.ndimage import gaussian_filter
    
    # Measure visit frequency
    visit_freq = np.zeros(raster_size)
    time_total = visit_freq.copy()
    for ii in list(range(len(walk_data['xinds']))):
        for jj in list(range(1, len(walk_data['xinds'][ii])-1)):
            # Running total of particles in this cell to find average
            visit_freq[walk_data['xinds'][ii][jj],
                       walk_data['yinds'][ii][jj]] += 1
            # Count the time in this cell as 0.5*last_dt + 0.5*next_dt
            last_dt = (walk_data['travel_times'][ii][jj] - \
                       walk_data['travel_times'][ii][jj-1])
            next_dt = (walk_data['travel_times'][ii][jj+1] - \
                       walk_data['travel_times'][ii][jj])
            time_total[walk_data['xinds'][ii][jj],
                       walk_data['yinds'][ii][jj]] += 0.5*(last_dt + next_dt)
    # Find mean time in each cell
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_time = time_total / visit_freq
    mean_time[visit_freq == 0] = 0
    # Prone to numerical outliers, so clip out extremes
    vmax = float(np.nanpercentile(mean_time, clip))
    mean_time = np.clip(mean_time, 0, vmax)
    
    # If applicable, do smoothing
    if sigma > 0:
        mean_time = gaussian_filter(mean_time, sigma=sigma)
        mean_time[mean_time==np.min(mean_time)] = np.nan
    else:
        mean_time[mean_time==0] = np.nan
    
    return mean_time


def unstruct2grid(coordinates,
                  quantity,
                  cellsize,
                  k_nearest_neighbors=3,
                  boundary=None,
                  crop=True):
    """Convert unstructured model outputs into gridded arrays.

    Interpolates model variables (e.g. depth, velocity) from an
    unstructured grid onto a Cartesian grid using inverse-distance-weighted
    interpolation. Assumes projected (i.e. "flat") geographic coordinates.
    Accepts coordinates in meters or decimal degrees. Extent of output
    rasters are based on extent of coordinates.
    (Function modeled after ANUGA plot_utils code)

    **Inputs** :

        coordinates : `list`
            List [] of (x,y) pairs or tuples of coordinates at which the
            interpolation quantities are located (e.g. centroids or vertices
            of an unstructured hydrodynamic model).

        quantity : `list`
            List [] of data to be interpolated with indices matching
            each (x,y) location given in coordinates. If quantity is
            depth, list would be formatted as [d1, d2, ... , dn].

        cellsize : `float or int`
            Length along one square cell face.

        k_nearest_neighbors : `int`, optional
            Number of nearest neighbors to use in the interpolation.
            If k>1, inverse-distance-weighted interpolation is used.

        boundary : `list`, optional
            List [] of (x,y) coordinates used to delineate the boundary of
            interpolation. Points outside the polygon will be assigned as nan.
            Format needs to match requirements of matplotlib.path.Path()

        crop : `bool`, optional
            If a boundary is specified, setting crop to True will eliminate
            any all-NaN borders from the interpolated rasters.

    **Outputs** :

        interp_func : `function`
            Nearest-neighbor interpolation function for gridding additional
            quantities. Quicker to use this output function on additional
            variables (e.g. later time-steps of an unsteady model) than
            to make additional function calls to unstruct2grid. Function
            assumes data have the same coordinates. It is used as follows:
            "new_gridded_quantity = interp_func(new_quantity)".

        gridded_quantity : `numpy.ndarray`
            Array of quantity after interpolation.

    """
    import matplotlib
    import scipy
    from scipy import interpolate

    cellsize = float(cellsize)

    # Make sure all input values are floats
    x = [float(i) for i, j in coordinates]
    y = [float(j) for i, j in coordinates]
    quantity = np.array([float(i) for i in quantity])
    if len(quantity) != len(x):
        raise ValueError("Coordinate and quantity arrays must be equal length")

    # Get some dimensions and make x,y grid
    nx = int(np.ceil((max(x)-min(x))/cellsize)+1)
    xvect = np.linspace(min(x), min(x)+cellsize*(nx-1), nx)
    ny = int(np.ceil((max(y)-min(y))/cellsize)+1)
    yvect = np.linspace(min(y), min(y)+cellsize*(ny-1), ny)

    gridX, gridY = np.meshgrid(xvect, yvect)

    inputXY = np.array([x[:], y[:]]).transpose()

    gridXY_array = np.array([np.concatenate(gridX),
                             np.concatenate(gridY)]).transpose()
    gridXY_array = np.ascontiguousarray(gridXY_array)
    
    # If a boundary has been specified, create array to index outside it
    if boundary is not None:
        path = matplotlib.path.Path(boundary)
        outside = ~path.contains_points(gridXY_array)

    # Create Interpolation function
    if k_nearest_neighbors == 1:  # Only use nearest neighbor
        index_qFun = interpolate.NearestNDInterpolator(inputXY,
                     np.arange(len(x), dtype='int64').transpose())
        gridqInd = index_qFun(gridXY_array)

        # Function to do the interpolation
        def interp_func(data):
            if isinstance(data, list):
                data = np.array(data)
            gridded_data = data[gridqInd].astype(float)
            if boundary is not None:
                gridded_data[outside] = np.nan # Crop to bounds
            gridded_data.shape = (len(yvect), len(xvect))
            gridded_data = np.flipud(gridded_data)
            if boundary is not None and crop is True:
                mask = ~np.isnan(gridded_data) # Delete all-nan border
                gridded_data = gridded_data[np.ix_(mask.any(1),
                                                   mask.any(0))]
            return gridded_data
    else:
        # Inverse-distance interpolation
        index_qFun = scipy.spatial.cKDTree(inputXY)
        NNInfo = index_qFun.query(gridXY_array, k=k_nearest_neighbors)
        # Weights for interpolation
        nn_wts = 1./(NNInfo[0]+1.0e-100)
        nn_inds = NNInfo[1]

        def interp_func(data):
            if isinstance(data, list):
                data = np.array(data)
            denom = 0.
            num = 0.
            for i in list(range(k_nearest_neighbors)):
                denom += nn_wts[:, i]
                num += data[nn_inds[:, i]].astype(float)*nn_wts[:, i]
            gridded_data = (num/denom)
            if boundary is not None:
                gridded_data[outside] = np.nan # Crop to bounds
            gridded_data.shape = (len(yvect), len(xvect))
            gridded_data = np.flipud(gridded_data)
            if boundary is not None and crop is True:
                mask = ~np.isnan(gridded_data) # Delete all-nan border
                gridded_data = gridded_data[np.ix_(mask.any(1),
                                                   mask.any(0))]
            return gridded_data

    # Finally, call the interpolation function to create array:
    gridded_quantity = interp_func(quantity)

    return interp_func, gridded_quantity
