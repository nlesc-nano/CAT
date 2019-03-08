Bond Dissociation Energy
========================

Calculate the bond dissociation energy (BDE) of ligands attached to the
surface of the core. The calculation consists of five distinct steps:

    1.  Dissociate all *n*2*(n-1)* combinations of 1 ligand (X), 1 Cd atom and
    1 other ligand (X).

    2.  Optimize the geometry of the CdX\ :sub:`2`\ structure with ADF MOPAC
    [1_, 2_, 3_].

    3.  Calculate the "electronic" contribution to the BDE (|dE|) with
    ADF MOPAC [1_, 2_, 3_] for all partially dissociated compounds
    created in step 1. This step consists of single point calculations.

    4.  Calculate the thermal contribution to the BDE (|ddG|) with
    ADF UFF [4_, 5_]. This step consists of geometry optimizations and
    frequency analyses.

    5.  |dG| = |dE_lvl1| + |ddG_lvl2| = |dE_lvl1| + ( |dG_lvl2| - |dE_lvl2| )



Default Settings
~~~~~~~~~~~~~~~~

::

    optional:
        qd:
            dissociate:
                job1: AMSJob
                job2: AMSJob
                s1: None
                s2: None

Arguments
~~~~~~~~~

**qd.dissociate.job1** |None|_, |str|_ or |type|_ = *None*

    The Job type used for calculating the "electronic" component (|dE_lvl1|) of
    the bond dissociation energy. Jobs c
    |None|_ will default to |AMSJob|_.

    |

**qd.dissociate.job2** |None|_, |bool|_, |str|_ or |type|_ = *None*

    The Job type used for calculating the thermal component (|ddG_lvl2|) of the
    bond dissociation energy. Setting it to |None|_ will default to |AMSJob|_.

    |

**qd.dissociate.s1** |None|_, |str|_ or |Settings|_ = *None*

    The Job type used for calculating the "electronic" component (|dE_lvl1|) of
    the bond dissociation energy.

    |

**qd.dissociate.s2** |None|_, |bool|_, |str|_ or |Settings|_ = *None*

    The Job type used for calculating the thermal component (|ddG_lvl2|) of the
    bond dissociation energy.

    |


.. _1: https://www.scm.com/doc/MOPAC/Introduction.html
.. _2: http://openmopac.net
.. _3: https://doi.org/10.1007/s00894-012-1667-x
.. _4: https://doi.org/10.1021/ja00051a040
.. _5: https://www.scm.com/doc/UFF/index.html

.. _AMSJob: https://www.scm.com/doc/plams/interfaces/ams.html#amsjob-api
.. _Job: https://www.scm.com/doc/plams/components/jobs.html#job-api
.. _Settings: https://www.scm.com/doc/plams/components/settings.html#api
.. _type: https://docs.python.org/3/library/functions.html#type
.. _bool: https://docs.python.org/3/library/stdtypes.html#boolean-values
.. _str: https://docs.python.org/3/library/stdtypes.html#str
.. _list: https://docs.python.org/3/library/stdtypes.html#list
.. _int: https://docs.python.org/3/library/functions.html#int
.. _None: https://docs.python.org/3/library/constants.html#None

.. |AMSJob| replace:: ``AMSJob``
.. |Job| replace:: ``Job``
.. |Settings| replace:: ``Settings``
.. |type| replace:: ``type``
.. |bool| replace:: ``bool``
.. |str| replace:: ``str``
.. |list| replace:: ``list``
.. |int| replace:: ``int``
.. |None| replace:: ``None``

.. |dE| replace:: d\ *E*
.. |dE_lvl1| replace:: d\ *E*\ :sub:`lvl1`
.. |dE_lvl2| replace:: d\ *E*\ :sub:`lvl2`
.. |dG| replace:: d\ *G*
.. |dG_lvl2| replace:: d\ *G*\ :sub:`lvl2`
.. |ddG| replace:: dd\ *G*
.. |ddG_lvl2| replace:: dd\ *G*\ :sub:`lvl2`
