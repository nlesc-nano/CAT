.. _Gettings Started:

General Overview & Getting Started
==================================

A basic recipe for running **CAT**:

1.  Create two directories named ‘core’ and ‘ligand’. The 'core' directory
should contain the input cores & the 'ligand' should contain the input
ligands. The quantum dots will be exported to the 'QD' directory.

2. 	Customize the job settings to your liking, see
CAT/examples/input_settings.yaml_ for an example.
Note: everything under the ``optional`` section does **not** have to be
included in the input settings.
As is implied by the name, everything in ``optional`` is completely optional.

3.  Run **CAT** with the following command:
``init_cat input_settings.yaml``

4.  Congratulations, you just ran
**CAT**!

The default **CAT** settings, at various levels of verbosity, are provided
below.

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O


Verbose default Settings
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O

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
            subset: null

        ligand:
            dirname: ligand
            optimize: True
            split: True
            functional_groups: null
            cosmo-rs: False

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            bulkiness: False
            activation_strain: False
            dissociate: False

Maximum verbose default Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O

    optional:
        database:
            dirname: database
            read: (core, ligand, qd)
            write: (core, ligand, qd)
            overwrite: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            dummy: Cl
            subset: null

        ligand:
            dirname: ligand
            split: True
            functional_groups: null
            cosmo-rs: False
            optimize:
                use_ff: False
                job1: null
                s1: null
                job2: null
                s2: null

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            bulkiness: False
            activation_strain: False
            dissociate:
                core_atom: Cd
                lig_count: 2
                keep_files: True
                core_core_dist: 5.0
                lig_core_dist: 5.0
                topology: {}

                job1: False
                s1: False
                job2: False
                s2: False

.. _input_settings.yaml: https://github.com/BvB93/CAT/blob/devel/examples/input_settings.yaml
