Further Optional Arguments
==========================

There are a number of arguments which can be used to modify the
functionality and behaviour of the quantum dot builder. Herein an
overview is provided.

Optional Arguments
~~~~~~~~~~~~~~~~~~

**dir_names** |list|_ [|str|_] = [*core*, *ligand*, *QD*]

    The names of the (to be created) folders.
    By default, ligand structures will be stored and read from dir_names[0],
    cores will be stored and read dir_names[1] and the combined core+ligands
    will be stored and read from dir_names[2]. Structures can be read from
    different folders if their filename is prepended with its absolute path,
    thus effectively ignoring this argument.

    |

**use_database** |bool|_ = *True*

    Enables or disables the storing and pulling of structures and properties
    from a user-created database (stored in .json and .xlsx formats).
    The script will attempt to pull a structure from the database if a match
    is found between a current input ligand and/or core+ligands and a previously
    optimized structure.

    |

**core.dummy** |str|_ or |int|_ = *Cl*

    The atomic number or atomic symbol of the atoms in the core that is to be
    replaced with ligands. Alternatively, dummy atoms can be manually specified
    with the core_indices variable.

    |

**ligand.optimize** |bool|_ = *True*

    Optimize the geometry of the to be attached ligands.
    The ligand is split into one or multiple (more or less) linear fragments, which
    are subsequently optimized (RDKit UFF [1_, 2_, 3_]) and reassembled while
    checking for the optimal dihedral angle. The ligand fragments are biased
    towards more linear conformations to minimize inter-ligand repulsion once the
    ligands are attached to the core.

    |

**ligand.split** |bool|_ = *True*

    If *False*: The ligand in its entirety is to be attached to the core.

        NR\ :sub:`4`\ :sup:`+` \                    -> NR\ :sub:`4`\ :sup:`+` \

        O\ :sub:`2`\CR                              -> O\ :sub:`2`\CR

        HO\ :sub:`2`\CR                             -> HO\ :sub:`2`\CR

        H\ :sub:`3`\CO\ :sub:`2`\CR                 -> H\ :sub:`3`\CO\ :sub:`2`\CR

    If *True*: A proton, counterion or functional group is to be removed from
    the ligand before attachment to the core.

        X\ :sup:`-`\.NR\ :sub:`4`\                  -> NR\ :sub:`4`\ :sup:`+` \

        HO\ :sub:`2`\CR                             -> O\ :sup:`-`\ :sub:`2`\CR

        Na\ :sup:`+`\.O\ :sup:`-`\ :sub:`2`\CR	    -> O\ :sup:`-`\ :sub:`2`\CR

        H\ :sub:`3`\CO\ :sub:`2`\CR                 -> O\ :sup:`-`\ :sub:`2`\CR

    |

**ligand.crs** |bool|_ = *False*

    Perform a property calculation with COSMO-RS [4_, 5_, 6_, 7_]; the COSMO
    surfaces are constructed using ADF MOPAC [8_, 9_, 10_].

    The following properties are calculated:

    1. The surface area of the ligand (A\ :sup:`2`\) as defined by its COSMO
    surface.

    2. The volume of the ligand (A\ :sup:`3`\) as defined by the volume
    encompassed by its COSMO surface.

    3. The solvation energy of the ligand (kcal mol\ :sup:`-1`\), at infinite
    dilution, in the following solvents: acetone, acetonitrile, dimethyl
    formamide (DMF), dimethyl sulfoxide (DMSO), ethyl acetate, ethanol,
    *n*-hexane, toluene and water.

    |

**qd.optimize** |bool|_ = *False*

    Optimize the quantum dot (i.e. core + all ligands) with ADF UFF [3_, 11_].
    The geometry of the core and ligand atoms directly attached to the core
    are frozen during this optimization.

    |

**qd.int** (|bool|_) = *False*

    Perform an activation strain analyses [12_, 13_, 14_] (kcal mol\ :sup:`-1`\)
    on the ligands attached to the quantum dot surface with RDKit UFF [1_, 2_, 3_].

    The core is removed during this process; the analyses is thus exclusively
    focused on ligand deformation and inter-ligand interaction.
    Yields three terms:

    1.  d\ *E*\ :sub:`strain`\  : 	The energy required to deform the ligand
    from their equilibrium geometry to the geometry they adopt on the quantum
    dot surface. This term is, by definition, destabilizing. Also known as the
    preperation energy (d\ *E*\ :sub:`prep`\).

    2.  d\ *E*\ :sub:`int`\  :	The mutual interaction between all deformed
    ligands. This term is characterized by the non-covalent interaction between
    ligands (UFF Lennard-Jones potential) and, depending on the inter-ligand
    distances, can be either stabilizing or destabilizing.

    3.  d\ *E* :	The sum of d\ *E*\ :sub:`strain`\  and d\ *E*\ :sub:`int`\ .
    Accounts for both the destabilizing ligand deformation and
    (de-)stabilizing interaction between all ligands in the absence of the core.

    |

**qd.bde** (|bool|_) = *False*


.. _1: http://www.rdkit.org
.. _2: https://github.com/rdkit/rdkit
.. _3: https://doi.org/10.1021/ja00051a040
.. _4: https://www.scm.com/doc/COSMO-RS/index.html
.. _5: https://doi.org/10.1021/j100007a062
.. _6: https://doi.org/10.1021/jp980017s
.. _7: https://doi.org/10.1139/V09-008
.. _8: https://www.scm.com/doc/MOPAC/Introduction.html
.. _9: http://openmopac.net
.. _10: https://doi.org/10.1007/s00894-012-1667-x
.. _11: https://www.scm.com/doc/UFF/index.html
.. _12: https://doi.org/10.1002/9780470125922.ch1
.. _13: https://doi.org/10.1002/wcms.1221
.. _14: https://doi.org/10.1021/acs.jpcc.5b02987

.. _bool: https://docs.python.org/3/library/stdtypes.html#boolean-values
.. _str: https://docs.python.org/3/library/stdtypes.html#str
.. _list: https://docs.python.org/3/library/stdtypes.html#list
.. _int: https://docs.python.org/3/library/functions.html#int

.. |bool| replace:: ``bool``
.. |str| replace:: ``str``
.. |list| replace:: ``list``
.. |int| replace:: ``int``
