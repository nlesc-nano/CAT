General Overview & Getting Started
==================================

A basic recipe for running **CAT**:**

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

Default Settings
~~~~~~~~~~~~~~~~

::

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O

More verbose default Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

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
            mol_format: [pdb, xyz]
            mongodb: False

        core:
            dirname: core
            dummy: Cl

        ligand:
            dirname: ligand
            optimize: True
            split: True
            cosmo-rs: False

        qd:
            dirname: QD
            optimize: False
            activation_strain: False
            dissociate: False

.. _input_settings.yaml: https://github.com/BvB93/CAT/blob/devel/examples/input_settings.yaml
