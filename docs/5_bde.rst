.. _Bond Dissociation Energy:

Bond Dissociation Energy
========================

Calculate the bond dissociation energy (BDE) of ligands attached to the
surface of the core. The calculation consists of five distinct steps:

    1.  Dissociate all combinations of *n* ligands (Y, see
    **qd.dissociate.lig_count**) and an atom from the core (X, see
    **qd.dissociate.core_atom**) within a radius *r* from aforementioned
    core atom (see **qd.dissociate.lig_core_dist** and
    **qd.dissociate.core_core_dist**).
    The dissociated compound has the general structure of |XYn|.

    2.  Optimize the geometry of |XYn| at the first level of theory
    (lvl1): ADF MOPAC [1_, 2_, 3_].

    3.  Calculate the "electronic" contribution to the BDE (|dE|)
    at the first level of theory (lvl1): ADF MOPAC [1_, 2_, 3_].
    This step consists of single point calculations of the complete
    quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

    4.  Calculate the thermalchemical contribution to the BDE (|ddG|) at the
    second level of theory (lvl2): ADF UFF [4_, 5_]. This step consists
    of geometry optimizations and frequency analyses of the same
    compounds used for step 3.

    5.  |dG| = |dE_lvl1| + |ddG_lvl2| = |dE_lvl1| + ( |dG_lvl2| - |dE_lvl2|
    ).

Default Settings
~~~~~~~~~~~~~~~~

::

    optional:
        qd:
            dissociate:
                core_atom: Cd
                lig_count: 2
                core_core_dist: 5.0
                lig_core_dist: 5.0
                topology:
                    7: vertice
                    8: edge
                    10: face

                job1: AMSJob
                s1: True
                job2: AMSJob
                s2: True

Arguments
~~~~~~~~~

**qd.dissociate.core_atom** |str|_ or |int|_ = *Cd*

    The atomic number or atomic symbol of the core atoms (X) which are to be
    dissociated. The core atoms are dissociated in combination with *n* ligands
    (Y, see **qd.dissociate.lig_count**).
    Yields a compound with the general formula |XYn|.

    |

**qd.dissociate.lig_count** |int|_ = *2*

    The number of ligands, *n*, which is to be dissociated in combination
    with a single core atom (X, see **qd.dissociate.core_atom**).
    Yields a compound with the general formula |XYn|.

    |

**qd.dissociate.core_core_dist** |float|_ = *5.0*

    The maximum to be considered distance (Ångström) between atoms in
    **qd.dissociate.core_atom**.
    Used for determining the topology of the core atom
    (see **qd.dissociate.topology**) and whether it is exposed to the
    surface of the core or not. It is recommended to use a radius which
    encapsulates a single (complete) shell of neighbours.

    |

**qd.dissociate.lig_core_dist** |float|_ = *5.0*

    Dissociate all possible combinations of *n* ligands and a single core atom
    (see **qd.dissociate.core_atom**) within a given radius (Ångström)
    from aforementioned core atom. The number of ligands dissociated in
    combination with a single core atom is controlled by
    **qd.dissociate.lig_count**.

    .. image:: _images/BDE_XY2.png
        :scale: 25 %
        :align: center

    |

**qd.dissociate.topology** |dict|_ =
{*7*: *vertice*, *8*: *edge*, *10*: *face*}

    A dictionary which translates the number neighbouring core atoms
    (see **qd.dissociate.core_atom** and **qd.dissociate.core_core_dist**)
    into a topology. Keys represent the number of neighbours, values represent
    the matching topology.

    Note: values can take on any user-specified value (*e.g.* Miller indices)
    and are thus not limited to *vertice*, *edge* and/or *face*.

    |

Arguments - Job Customization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**qd.dissociate.job1** |type|_, |str|_ or |bool|_ = *AMSJob*

    A |type|_ object of a |Job|_ subclass, used for calculating the
    "electronic" component (|dE_lvl1|) of the bond dissociation energy.
    Involves single point calculations.

    Alternatively, an alias (|str|_) can be provided for a specific
    job type (see :ref:`Type Aliases`).

    Setting it to *True* (|bool|_) will default to |type|_ (|AMSJob|_),
    while *False* (|bool|_) is equivalent to
    ``optional.qd.dissociate = False``.

    |

**qd.dissociate.s1** |Settings|_, |str|_ or |bool|_ =

    ::

        s1:
            input:
                mopac:
                    model: PM7
                ams:
                    system:
                        charge: 0

    The job |Settings|_ used for calculating the "electronic" component
    (|dE_lvl1|) of the bond dissociation energy.

    Alternatively, a path (|str|_) can be provided to .json or .yaml file
    containing the job settings.

    Setting it to *True* (|bool|_) will default to the *MOPAC* block in
    CAT/data/templates/qd.yaml_, while *False* (|bool|_) is equivalent to
    ``optional.qd.dissociate = False``.

    |

**qd.dissociate.job2** |type|_, |str|_ or |bool|_ = *AMSJob*

    A |type|_ object of a |Job|_ subclass, used for calculating the
    thermal component (|ddG_lvl2|) of the bond dissociation energy.
    Involves a geometry reoptimizations and frequency analyses.

    Alternatively, an alias (|str|_) can be provided for a specific
    job type (see :ref:`Type Aliases`).


    Setting it to *True* (|bool|_) will default to |type|_ (|AMSJob|_),
    while *False* (|bool|_) will skip the thermochemical analysis completely.

    |

**qd.dissociate.s2** |Settings|_, |str|_ or |bool|_ =

    ::

        s2:
            input:
                uff:
                    library: uff
                ams:
                    system:
                        charge: 0
                        bondorders:
                            _1: null

    The job |Settings|_ used for calculating the thermal component (|ddG_lvl2|)
    of the bond dissociation energy.

    Alternatively, a path (|str|_) can be provided to .json or .yaml file
    containing the job settings.

    Setting it to *True* (|bool|_) will default to the the *MOPAC* block in
    CAT/data/templates/qd.yaml_, while *False* (|bool|_) will skip the
    thermochemical analysis completely.

    |

.. _1: https://www.scm.com/doc/MOPAC/Introduction.html
.. _2: http://openmopac.net
.. _3: https://doi.org/10.1007/s00894-012-1667-x
.. _4: https://doi.org/10.1021/ja00051a040
.. _5: https://www.scm.com/doc/UFF/index.html
.. _qd.yaml: https://github.com/BvB93/CAT/blob/master/CAT/data/templates/qd.yaml

.. _AMSJob: https://www.scm.com/doc/plams/interfaces/ams.html#amsjob-api
.. _Job: https://www.scm.com/doc/plams/components/jobs.html#job-api
.. _Settings: https://www.scm.com/doc/plams/components/settings.html#api
.. _type: https://docs.python.org/3/library/functions.html#type
.. _bool: https://docs.python.org/3/library/stdtypes.html#boolean-values
.. _str: https://docs.python.org/3/library/stdtypes.html#str
.. _list: https://docs.python.org/3/library/stdtypes.html#list
.. _dict: https://docs.python.org/3/library/stdtypes.html#dict
.. _int: https://docs.python.org/3/library/functions.html#int
.. _float: https://docs.python.org/3/library/functions.html#float
.. _None: https://docs.python.org/3/library/constants.html#None

.. |AMSJob| replace:: ``AMSJob``
.. |Job| replace:: ``Job``
.. |Settings| replace:: ``Settings``
.. |type| replace:: ``type``
.. |bool| replace:: ``bool``
.. |str| replace:: ``str``
.. |list| replace:: ``list``
.. |dict| replace:: ``dict``
.. |int| replace:: ``int``
.. |float| replace:: ``float``
.. |None| replace:: ``None``

.. |dE| replace:: d\ *E*
.. |dE_lvl1| replace:: d\ *E*\ :sub:`lvl1`
.. |dE_lvl2| replace:: d\ *E*\ :sub:`lvl2`
.. |dG| replace:: d\ *G*
.. |dG_lvl2| replace:: d\ *G*\ :sub:`lvl2`
.. |ddG| replace:: dd\ *G*
.. |ddG_lvl2| replace:: dd\ *G*\ :sub:`lvl2`
.. |XYn| replace:: XY\ :sub:`n`
.. |Yn| replace:: Y\ :sub:`n`
