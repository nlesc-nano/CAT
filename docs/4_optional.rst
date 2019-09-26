.. _Optional:

Optional
========

There are a number of arguments which can be used to modify the
functionality and behaviour of the quantum dot builder. Herein an
overview is provided.

Note: Inclusion of this section in the input file is not required,
assuming one is content with the default settings.

Index
~~~~~

========================================= ==================================================================================
Option                                    Description
========================================= ==================================================================================
:attr:`optional.database.dirname`         The name of the directory where the database will be stored.
:attr:`optional.database.read`            Attempt to read results from the database before starting calculations.
:attr:`optional.database.write`           Export results to the database.
:attr:`optional.database.overwrite`       Allow previous results in the database to be overwritten.
:attr:`optional.database.mol_format`      The file format(s) for exporting moleculair structures.
:attr:`optional.database.mongodb`         Options related to the MongoDB format.

:attr:`optional.core.dirname`             The name of the directory where all cores will be stored.
:attr:`optional.core.dummy`               Atomic number of symbol of the core dummy atoms.

:attr:`optional.ligand.dirname`           The name of the directory where all ligands will be stored.
:attr:`optional.ligand.optimize`          Optimize the geometry of the to-be attached ligands.
:attr:`optional.ligand.functional_groups` Manually specify SMILES strings representing functional groups.
:attr:`optional.ligand.split`             If the ligand should be attached in its entirety to the core or not.
:attr:`optional.ligand.cosmo-rs`          Perform a property calculation with COSMO-RS on the ligand.
:attr:`optional.ligand.bulkiness`         Perform a ligand bulkiness calculation.

:attr:`optional.qd.dirname`               The name of the directory where all quantum dots will be stored.
:attr:`optional.qd.optimize`              Optimize the quantum dot (i.e. core + all ligands) .
:attr:`optional.qd.activation_strain`     Perform an activation strain analyses.
:attr:`optional.qd.dissociate`            Calculate the ligand dissociation energy.
========================================= ==================================================================================

Default Settings
~~~~~~~~~~~~~~~~

.. code::

    optional:
        database:
            dirname: database
            read: True
            write: True
            overwrite: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            dummy: Cl

        ligand:
            dirname: ligand
            optimize: True
            functional_groups: null
            split: True
            cosmo-rs: False
            bulkiness: False

        qd:
            dirname: qd
            optimize: False
            activation_strain: False
            dissociate: False


Arguments
~~~~~~~~~

Database
--------

.. attribute:: optional.database

    All database-related settings.

    .. note::
        For :attr:`optional.database` settings to take effect the `Data-CAT <https://github.com/nlesc-nano/data-CAT>`_ package has to be installed.

    Example:

    .. code::

        optional:
            database:
                dirname: database
                read: True
                write: True
                overwrite: False
                mol_format: (pdb, xyz)
                mongodb: False

|

    .. attribute:: optional.database.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default Value** - ``"database"``

        The name of the directory where the database will be stored.

        The database directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.database.read

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Attempt to read results from the database before starting calculations.

        Before optimizing a structure, check if a geometry is available from
        previous calculations. If a match is found, use that structure and
        avoid a geometry reoptimizations. If one wants more control then the
        boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures will be read only for a
        specific subset.


        .. admonition:: Example

            Example #1:

            .. code::

                optional:
                    database:
                        read: (core, ligand, qd)  # This is equivalent to read: True

            Example #2:

            .. code::

                optional:
                    database:
                        read: ligand


    .. attribute:: optional.database.write

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Export results to the database.

        Previous results will **not** be overwritten unless
        :attr:`optional.database.overwrite` = ``True``. If one wants more control then
        the boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures written for for a specific
        subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.overwrite

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``False``

        Allow previous results in the database to be overwritten.

        Only apllicable if :attr:`optional.database.write` = ``True``.
        If one wants more control then the boolean can be substituted for
        a list of strings (*i.e.* ``"core"``, ``"ligand"`` and/or ``"qd"``), meaning
        that structures written for for a specific subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.mol_format

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("pdb", "xyz")``

        The file format(s) for exporting moleculair structures.

        By default all structures are stored in the .hdf5 format as
        (partially) de-serialized .pdb files. Additional formats can be
        requisted with this keyword.
        Accepted values: ``"pdb"``, ``"xyz"``, ``"mol"`` and/or ``"mol2"``.


    .. attribute:: optional.database.mongodb

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default Value** – ``False``

        Options related to the MongoDB format.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Database`:.

|

Core
----

.. attribute:: optional.core

    All settings related to the core.

    Example:

    .. code::

        optional:
            core:
                dirname: core
                dummy: Cl

|

    .. attribute:: optional.core.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"core"``

        The name of the directory where all cores will be stored.

        The core directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.core.dummy

        :Parameter:     * **Type** - :class:`str` or :class:`int`
                        * **Default value** – ``17``


        Atomic number of symbol of the core dummy atoms.

        The atomic number or atomic symbol of the atoms in the core which are to be
        replaced with ligands. Alternatively, dummy atoms can be manually specified
        with the core_indices variable.

|

Ligand
------

.. attribute:: optional.ligand

    All settings related to the ligands.

    Example:

    .. code::

        optional:
            ligand:
                dirname: ligand
                optimize: True
                functional_groups: null
                split: True
                cosmo-rs: False
                bulkiness: False

|

    .. attribute:: optional.ligand.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"ligand"``

        The name of the directory where all ligands will be stored.

        The ligand directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.ligand.optimize

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Optimize the geometry of the to-be attached ligands.

        The ligand is split into one or multiple (more or less) linear fragments,
        which are subsequently optimized (RDKit UFF [1_, 2_, 3_]) and reassembled
        while checking for the optimal dihedral angle. The ligand fragments are
        biased towards more linear conformations to minimize inter-ligand
        repulsion once the ligands are attached to the core.


    .. attribute:: optional.ligand.functional_groups

        :Parameter:     * **Type** - :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** – ``None``

        Manually specify SMILES strings representing functional groups.

        For example, with :attr:`optional.ligand.functional_groups` = ``("O[H]", "[N+].[Cl-]")`` all
        ligands will be searched for the presence of hydroxides and ammonium chlorides.

        The first atom in each SMILES string (*i.e.* the "anchor") will be used for attaching the ligand
        to the core, while the last atom (assuming :attr:`optional.ligand.split` = ``True``) will be
        dissociated from the ligand and disgarded.

        If not specified, the default functional groups of **CAT** are used.

        .. note::
            This argument has no value be default and will thus default to SMILES strings of the default
            functional groups supported by **CAT**.

        .. note::
            The yaml format uses ``null`` rather than ``None`` as in Python.

    .. attribute:: optional.ligand.split

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        If ``False``: The ligand is to be attached to the core in its entirety .

        =================== ==================
        Before              After
        =================== ==================
        :math:`{NR_4}^+`    :math:`{NR_4}^+`
        :math:`O_2 CR`      :math:`O_2 CR`
        :math:`HO_2 CR`     :math:`HO_2 CR`
        :math:`H_3 CO_2 CR` :math:`H_3 CO_2 CR`
        =================== ==================

        ``True``: A proton, counterion or functional group is to be removed from
        the ligand before attachment to the core.

        ========================= ==================
        Before                    After
        ========================= ==================
        :math:`Cl^- + {NR_4}^+`   :math:`{NR_4}^+`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`Na^+ + {O_2 CR}^-` :math:`{O_2 CR}^-`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`H_3 CO_2 CR`       :math:`{O_2 CR}^-`
        ========================= ==================


    .. attribute:: optional.ligand.cosmo-rs

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``


        Perform a property calculation with COSMO-RS [4_, 5_, 6_, 7_] on the ligand.

        The COSMO surfaces are by default constructed using ADF MOPAC [8_, 9_, 10_].

        The solvation energy of the ligand and its activity coefficient are
        calculated in the following solvents: acetone, acetonitrile,
        dimethyl formamide (DMF), dimethyl sulfoxide (DMSO), ethyl acetate,
        ethanol, *n*-hexane, toluene and water.


    .. attribute:: optional.ligand.bulkiness

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``


        Calculate the ligand dissociation energy.

        Given a set of angles :math:`\phi`, the bulkiness factor :math:`V_{bulk}` is defined below.
        Angles are constructed according to :math:`\phi = \angle_{ABC}`,
        where :math:`A` represents a set of all ligand atoms,
        :math:`B` is the ligand anchor atom and
        :math:`C` is the mean position of all ligand atoms (*i.e.* the ligand center).

        .. math::
            W_{bulk} = \frac{1}{n} \sum_{i}^{n} e^{\sin \phi_{i}}

        Conceptually, the bulkiness factor :math:`V_{bulk}` is related to ligand (half-)cones angles,
        with the key difference that :math:`V_{bulk}` builds on top of it,
        representing an estimate of mean inter-ligand steric interactions.

        See also
        --------
        `Ligand cone angle <https://en.wikipedia.org/wiki/Ligand_cone_angle>`_:
            The ligand cone angle is a measure of the steric bulk of a ligand in
            a transition metal complex.

|

QD
--

.. attribute:: optional.qd

    All settings related to the quantum dots.

    Example:

    .. code::

        optional:
            qd:
                dirname: QD
                optimize: False
                activation_strain: False
                dissociate: False

|

    .. attribute:: optional.qd.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"qd"``

        The name of the directory where all quantum dots will be stored.

        The quantum dot directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.qd.optimize

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Optimize the quantum dot (i.e. core + all ligands) .

        By default the calculation is performed with ADF UFF [3_, 11_].
        The geometry of the core and ligand atoms directly attached to the core
        are frozen during this optimization.


    .. attribute:: optional.qd.activation_strain

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Perform an activation strain analyses [12_, 13_, 14_].

        The activation strain analyses (kcal mol\ :sup:`-1`\) is performed
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
        Accounts for both the destabilizing ligand deformation and (de-)stabilizing
        interaction between all ligands in the absence of the core.


    .. attribute:: optional.qd.dissociate

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Calculate the ligand dissociation energy.

        Calculate the ligand dissociation energy (BDE) of ligands attached to the
        surface of the core. See :ref:`Bond Dissociation Energy` for more details.
        The calculation consists of five distinct steps:

            1.  Dissociate all combinations of |n| ligands (|Y|) and an atom from the core (|X|)
            within a radius *r* from aforementioned core atom.
            The dissociated compound has the general structure of |XYn|.

            2.  Optimize the geometry of |XYn| at the first level of theory
            (:math:`1`). Default: ADF MOPAC [1_, 2_, 3_].

            3.  Calculate the "electronic" contribution to the BDE (|dE|)
            at the first level of theory (:math:`1`): ADF MOPAC [1_, 2_, 3_].
            This step consists of single point calculations of the complete
            quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

            4.  Calculate the thermalchemical contribution to the BDE (|ddG|) at the
            second level of theory (:math:`2`). Default: ADF UFF [4_, 5_]. This step
            consists of geometry optimizations and frequency analyses of the same
            compounds used for step 3.

            5.  :math:`\Delta G_{tot} = \Delta E_{1} + \Delta \Delta G_{2} = \Delta E_{1} + (\Delta G_{2} - \Delta E_{2})`.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Bond Dissociation Energy`:.



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

.. |dE| replace:: :math:`\Delta E`
.. |dE_lvl1| replace:: :math:`\Delta E_{1}`
.. |dE_lvl2| replace:: :math:`\Delta E_{2}`
.. |dG| replace:: :math:`\Delta G_{tot}`
.. |dG_lvl2| replace:: :math:`\Delta G_{2}`
.. |ddG| replace:: :math:`\Delta \Delta G`
.. |ddG_lvl2| replace:: :math:`\Delta \Delta G_{2}`
.. |XYn| replace:: :math:`XY_{n}`
.. |Yn| replace:: :math:`Y_{n}`
.. |n| replace:: :math:`{n}`
.. |X| replace:: :math:`X`
.. |Y| replace:: :math:`Y`
