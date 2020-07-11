.. _background:

==========
Background
==========

Motivation
----------
particlerouting was developed to provide an open-source and general method for routing passive tracer particles in a landscape scale flow field. Methods to simulate Lagrangian particle transport exist in adjacent fields and have been developed to solve problems at large scales (OpenDrift [1]) and at small scales (freshkiss3d [2]). The popular Delft3D suite of modules contains a method for Lagrangian particle tracking called Delft3D-PART [3], however as part of the larger, complex Delft3D ecosystem, this code is harder to interrogate and understand directly. 

By developing particlerouting as a standalone python package, we hope to simplify the simulation of Lagrangian particle transport so that it may be applied to the results from any hydrodynamic simulation. 

Theory
------
Particle transport is simulated through the use of a weighted random walk framework [4]. In particular, we apply methods used in the numerical model DeltaRCM [5,6], to weight the random walk by the local water slopes and the water discharge components.

particlerouting assumes a rectangular grid and routes particles in a D-8 fashion assuming that the potential future location of a given particle is one of the surrounding 8 cells. 

The weights for the random walk are sensitive to 2 parameters, :math:`{\gamma}` and :math:`{\theta}`.

The routing direction F* is the estimate of the local downstream of the flow [5]. This direction F* is comprised of :math:`{F_{sfc}}`, and :math:`{F_{int}}`, calculated based on the water surface gradient and the discharge respectively. The proportional combination of these directional components is dictated by the parameter :math:`{\gamma}`:

.. math::
   :nowrap:

   \begin{eqnarray}
      F^{*} = \gamma F_{sfc} + (1-\gamma) F_{int}
   \end{eqnarray}

In this way, the :math:`{\gamma}` parameter controls the proportional dependence of water surface slope and water discharge on the downstream direction. A :math:`{\gamma}` of 0 means that the water slope is ignored when determining the downstream flow direction, and a :math:`{\gamma}`: of 1 means that the discharge values (the flow field) is ignored.

The second weighting parameter, :math:`{\theta}` modifies the routing weight for the random walk based on the local water depth value. The routing weights for each neighboring cell are calculated per [5], using the following equation:

.. math::
   :nowrap:

   \begin{eqnarray}
      w_i = \frac{\frac{1}{R_i} \text{max} \left(0, F \cdot d_i \right)}{\Delta_i}
   \end{eqnarray}

The resistance value, :math:`{R_i}`, is computed using the local water depth, :math:`{h_i}` and the parameter :math:`{\theta}`:

.. math::
   :nowrap:

   \begin{eqnarray}
      R_i = \frac{1}{h_i^\theta}
   \end{eqnarray}

When :math:`{\theta}` is 0, then the local depth values are not impacting the weighting factor associated with each cell. As :math:`{\theta}` gets larger, the local routing weights have an increasingly large dependence on the water depth of the neighboring cells.

References
----------
[1] Dagestad, K.-F., Röhrs, J., Breivik, Ø., and Ådlandsvik, B.: OpenDrift v1.0: a generic framework for trajectory modelling, Geosci. Model Dev., 11, 1405-1420, https://doi.org/10.5194/gmd-11-1405-2018, 2018.

[2] ANGE Team, Freshkiss3D home page. Available at: https://freshkiss3d.gforge.inria.fr/

[3] Hydraulics, D. (2007). Delft3D-PART user manual version 2.13. WL| Delft Hydraulics, Delft.

[4] Pearson, K. (1905). The problem of the random walk. Nature, 72(1867), 342-342.

[5] Liang, M., Voller, V. R., & Paola, C. (2015). A reduced-complexity model for river delta formation-Part 1: Modeling deltas with channel dynamics. Earth Surface Dynamics, 3(1).

[6] Liang, M., Geleynse, N., Edmonds, D. A., & Passalacqua, P. (2015). A reduced-complexity model for river delta formation-Part 2: Assessment of the flow routing scheme. Earth Surface Dynamics, 3(1), 87.


