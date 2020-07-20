.. Quick Demo 1_:

Demo 1- Using the High-Level API
--------------------------------
In this demo, we show how particles can be routed along a simulated river delta, example output is from the delta simulation model, `DeltaRCM <https://github.com/DeltaRCM/pyDeltaRCM_WMT>`_. To do this, we will be using one of the high-level API functions provided in :download:`routines.py <../../../particlerouting/routines.py>`.

First we load the sample parameters.

.. doctest::

    >>> from particlerouting.example_data import define_params as dp
    >>> import particlerouting as pr
    >>> rcmparams = dp.make_rcm_params()

We can visualize the water depth for this scenario from the parameters. If you'd like to download this portion of the demo as a standalone script, it is available :download:`here <../pyplots/quickstart/demo1_depth.py>`.

.. doctest::

    >>> plt.imshow(rcmparams.depth)
    >>> plt.colorbar()
    >>> plt.title('Water Depth')
    >>> plt.show()

.. plot:: quickstart/demo1_depth.py

Now let's route 50 particles for 50 iterations.

.. doctest::

    >>> pr.routines.steady_plots(rcmparams, 50, 'demo-1')

The `steady_plots()` function saves plots of each iteration of the particle movement to the subfolder 'demo-1/figs'. The final plot is shown below; the initial particle locations are shown as blue dots, and the final particle locations are red dots. If you'd like to download this portion of the demo as a standalone script, it is available :download:`here <../pyplots/quickstart/demo1.py>`.

.. plot:: quickstart/demo1.py
