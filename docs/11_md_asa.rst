Ensemble-Averaged Activation Strain Analysis
============================================
.. math::
    :label: 1a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}

Herein we describe an Ensemble-Averaged extension of the
activation/strain  analysis (ASA; also known as the
distortion/interaction model), wherein the ASA is utilized
for the analyses of entire molecular dynamics trajectories.
The implementation utilizes CHARMM-style forcefields for the
calculation of all energy terms.

.. note::
    Throughout this document an overline will be used to distinguish between "normal"
    and ensemble-averaged quantities: *e.g.* :math:`E_{\text{strain}}` versus
    :math:`\overline{E}_{\text{strain}}`.

|

Strain/Distortion
-----------------
The ensemble averaged strain :math:`\Delta \overline{E}_{\text{strain}}`
represents the distortion of all ligands with respect to their equilibrium
geometry.
Given an MD trajectory with :math:`m` iterations and :math:`n` ligands per
quantum dot, the energy is averaged over all :math:`m` MD iterations and
summed over all :math:`n` ligands.

The magnitude of this term is determined by all covalent and non-covalent
intra-ligand interactions.
As this term quantifies the deviation of a ligand from its equilibrium geometry,
it is, by definition, always positive.

.. math::
    :label: 2a

    \Delta E_{\text{strain}} = E_{\text{lig-pert}} - E_{\text{lig-eq}}
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{strain}} = \frac{1}{m} \sum_{i=0}^{m} \sum_{j=0}^{n}
    E_{\text{lig-pert}}(i, j) - E_{\text{lig-eq}}


.. math::
    :label: 3a

    \Delta E_{\text{strain}} = \Delta V_{\text{bond}} + \Delta V_{\text{angle}} +
    \Delta V_{\text{Urey-Bradley}} + \Delta V_{\text{dihedral}} + \Delta V_{\text{improper}} +
    \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`E_{\text{lig-eq}}` is herein the total energy of a (single) ligand at
its equilibrium geometry, while :math:`E_{\text{lig-pert}}(i, j)` is the
total energy of the (perturbed) ligand :math:`j` at MD iteration :math:`i`.

|

Interaction
-----------
The ensemble averaged interaction :math:`\Delta \overline{E}_{\text{int}}`
represents the mutual interaction between all ligands in a molecule.
The interaction is, again, averaged over all MD iterations and
summed over all ligand-pairs.

The magnitude of this term is determined by all non-covalent inter-ligand
interactions and can be either positive (dominated by Pauli and/or
Coulombic repulsion) or negative (dominated by dispersion and/or Coulombic
attraction).

.. math::
    :label: 4a

    \Delta E_{\text{int}} = \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (j, k)
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{int}} = \frac{1}{m} \sum_{i=0}^{m}
    \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (i, j, k)

.. math::
    :label: 5a

    \Delta E_{\text{int}} = \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`\Delta E_{\text{lig-int}}(i, j, k)` represents the pair-wise
interactions between ligands :math:`j` and :math:`k` at MD iteration :math:`i`.
Double counting is avoided by ensuring that :math:`k > j`.

.. note::
    In order to avoid the substantial Coulombic repulsion between negatively charged ligands,
    its parameters are substituted with those from its neutral (*i.e.* protonated) counterpart.
    This correction is applied, exclusively, for the calculation of :math:`\Delta E_{\text{lig-int}}`.

|

Total Energy
------------
The total (ensemble-averaged) energy is the sum of
:math:`\Delta \overline{E}_{\text{strain}}` and
:math:`\Delta \overline{E}_{\text{int}}`.
Note that the energy is associated with a set of :math:`n` ligands,
*i.e.* the distortion and mutual interaction between all :math:`n` ligands.
Division by :math:`n` will thus yield the averaged energy per ligand
per MD iteration.

.. math::
    :label: 6a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}
    = \frac{1}{m} \sum_{i=0}^{m} \Delta E_{\text{strain}}(i) + \Delta E_{\text{int}}(i)
