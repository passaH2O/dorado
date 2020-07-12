.. _example06:

Example 6 - Set Travel Time Target
==================================

In this example we will revisit the `ANUGA <https://github.com/GeoscienceAustralia/anuga_core>`_ model output we looked at in :ref:`example01`. This time, however, we will set a particle travel time target as opposed to a number of iterations for the particles to move.

In this method, we allow the particles to take as many iterations as they need (individually) to get as close as possible to the target travel time that we prescribe.

Full example script available :download:`here <../../../examples/set_timestep_anuga_particles.py>`.

After loading the data and establishing the parameters in the same way we did in :ref:`example01`, we will change the grid size so that the range of travel times we can obtain in this small example is a bit better constrained.

.. doctest::

   >>> params.dx = 10.

We note that by discretizing our grid and the particle routing procedure, we are unable to prescribe precise particle travel times. Due to the gridded scheme, we can only measure particle travel times when particles are located at the center of a cell. So when a target travel time is set, the particle travel times are not precisely that number, but rather the closest times corresponding to when each particle is in a cell center.

To set a travel time target, we will use the lower-level API and access the `particle_track` module directly, and then use some routines to plot the initial and final particle locations. The target travel time we will be using is 2100 seconds.

.. doctest::

   >>> walk_data = particle.run_iteration(target_time=2100)
   >>> pr.routines.plot_initial(anugaparams.depth, walk_data)
   >>> pr.routines.plot_final(anugaparams.depth, walk_depth)
   >>> plt.title('Initial and Final Particle Locations')
   >>> plt.show()

.. plot:: quickstart/demo2.py

To see how close we are to the prescribed target travel time of 2100 seconds, we will loop through the particles and get a list of the final travel times that we then print to the console.

.. doctest::

   >>> finaltimes = []
   >>> for i in list(range(0, params.Np_tracer)):
   >>>    finaltimes.append(walk_data['travel_times'][i][-1])
   >>> print('List of particle travel times for final particle locations: ' +
   >>>       str(np.round(finaltimes)))

   List of particle travel times for final particle locations: [2158.
   2118. 2088. 2149. 2096. 2206. 2104. 2104. 2207. 2137. 2097. 2147.
   2032. 2118. 2084. 2037. 2066. 2076. 2104. 2108. 2033. 2101. 2080.
   2072. 2031. 2041. 2076. 2063. 2125. 2102. 2140. 2178. 2173. 2097.
   2104. 2189. 2061. 2112. 2074. 2095. 2100. 2177. 2069. 2032. 2050.
   2086. 2036. 2109. 2078. 2047.]
