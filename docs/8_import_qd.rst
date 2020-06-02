.. _import_qd:

Importing Quantum Dots
======================

*WiP*: Import pre-built quantum dots rather than constructing them from scratch.


Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    input_qd:
        - Cd68Se55_ethoxide.xyz:
            ligand_smiles: '[O-]CC'
            ligand_anchor: '[O-]'

|

Arguments
~~~~~~~~~

.. attribute:: ligand_smiles

    :Parameter:     * **Type** - :class:`str`
                    * **Default value** – ``None``

    A SMILES string representing the ligand.
    The provided SMILES string will be used for identifying the core and all ligands.

    .. warning::
        This argument has no value be default and thus *must* be provided by the user.


.. attribute:: ligand_anchor

    :Parameter:     * **Type** - :class:`str`
                    * **Default value** – ``None``

    A SMILES string representing the achor functional group of the ligand.
    If the provided SMILES string consists of multiple atoms
    (*e.g.* a carboxylate: ``"[O-]C=O"``), than the first atom will be treated as anchor (``"[O-]"``).

    .. warning::
        This argument has no value be default and thus *must* be provided by the user.
