.. _example10:

Example 10 - Routing Particles in Parallel
==========================================

.. note:: This functionality has only been tested on UNIX-based systems.

In this example, we will use the parallel routing functionality to distribute the particle routing procedure across multiple CPU cores. Then we will compare the time it takes to route particles in serial against the time it takes to route an equivalent number of particles on 2 CPU cores. To do this, we use the native `multiprocessing` library that is provided with Python.

Full example script available :download:`here <../../../examples/parallel_comparison.py>`.

We will be running the parallel comparison on the ANUGA example domain used in :ref:`example01`. The only difference is that we will use a larger number of particles so that the impact of the parallelization can be measured even on this relatively small and simple domain.

So after initializing the parameters as we did in :ref:`example01`, we will define the number of particles and where they should be seeded.

.. doctest::

   >>> Np_tracer = 200
   >>> seed_xloc = list(range(20, 30))
   >>> seed_yloc = list(range(48, 53))

Next we will do the parallel routing using 2 CPU cores using the `parallel_routing` function. In this case, we are routing 200 particles on each core, and iterating each particle 50 times.

.. doctest::

   >>> particles = Particles(params)
   >>> from dorado.parallel_routing import parallel_routing
   >>> par_result = parallel_routing(particles, 50, Np_tracer, seed_xloc,
                                     seed_yloc, 2)

For comparison, we will run the comparable serial process twice so that the same number of particles are simulated.

.. doctest::

   >>> for z in list(range(0, 2)):
   >>>    particle = Particles(params)
   >>>    particle.generate_particles(Np_tracer, seed_xloc, seed_yloc)
   >>>    for i in list(range(0, 50)):
   >>>        all_walk = particle.run_iteration()

If we were to time this example (timing included in full example script), a roughly 2x speed-up can be observed. Speed-up may vary based on a variety of factors.

.. doctest::

   Serial Compute Time: 3.2187864780426025
   Parallel Compute Time: 1.6109538078308105
