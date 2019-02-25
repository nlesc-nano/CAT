General Overview & Getting Started
==================================

A basic recipe for running **CAT**:**

1.  Create two directories named ‘core’ and ‘ligand’. The 'core' directory should contain the input cores & the 'ligand' should contain the input ligands. The quantum dots will be exported to the 'QD' directory.

2. 	Customize the job settings to your liking, see CAT/examples/input_settings.yaml_ for an example.

3.  Run **CAT** with ``init_cat input_settings.yaml``

4.  Congratulations, you just ran **CAT**!

.. _qd-example: https://github.com/SCM-NV/qmflows/blob/master/test/QD_input_examples

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

    optional:
        dir_names: [core, ligand, QD]
        use_database: False
        core:
            dummy: Cl
        ligand:
            optimize: True
            split: True
            cosmo-rs: False
        qd:
            optimize: False
            activation_strain: False
            dissociate: False

.. _input_settings.yaml: https://github.com/BvB93/CAT/blob/devel/examples/input_settings.yaml