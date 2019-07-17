.. _Input Cores and Ligands:

input_cores & input_ligands
===========================

Thia section related relates the importing and processing of cores and ligands.
Ligand & cores can be imported from a wide range of different files and files
types, which can roughly be divided into three categories:

1.  Files containing coordinates of a single molecule: .xyz, .pdb & .mol files.
2.  Python objects: :class:`plams.Molecule`, :class:`rdkit.Chem.Mol` & SMILES strings (:class:`str`).
3.  Containers with one or multiple input molecules: directories & .txt files.

In the later case, the container can consist of multiple SMILES strings or
paths to .xyz, .pdb and/or .mol files. If necessary, containers are searched
recursively. Both absolute and relative paths are explored.

Default Settings
~~~~~~~~~~~~~~~~

.. code::

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O
        - OC(CCC)=O
        - OC(CCCC)=O

Optional arguments
~~~~~~~~~~~~~~~~~~

.. attribute:: .guess_bonds

    :Parameter:     * **Type** - :class:`bool`
                    * **Default value** – ``False``

    Try to guess bonds and bond orders in a molecule based on the types atoms
    and the relative of atoms. Is set to False by default, with the exception
    of .xyz files.

.. attribute:: .column

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The column containing the to be imported molecules.
    Relevant when importing structures from .txt and .xlsx files with
    multiple columns. Numbering starts from 0.


.. attribute:: .row

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The first row in a column which contains a molecule.
    Useful for when, for example, the very first row contains the title of
    aforementioned row, in which case row = 1 would be a sensible choice.
    Relevant for .txt and .xlsx files. Numbering starts from 0.

.. attribute:: .indices

    :Parameter:     * **Type** - :class:`int` or :class:`tuple` [:class:`int`]
                    * **Default value** – ``None``

    input_cores:
        Manually specify the atomic index of one ore more atom(s) in the core that
        will be replaced with ligands. If left empty, all atoms of a user-specified
        element (see :attr:`optional.cores.dummy`) will be replaced with
        ligands.

    input_ligands:
        Manually specify the atomic index of the ligand atom that will be attached
        to core (implying argument_dict: ``optional.ligand.split = False``).If two
        atomic indices are rovided, the bond between :class:`tuple` [``0``] and
        :class:`tuple` [``1``] will be broken and the molecule containing
        :class:`tuple` [``0``] is attached to the core,
        (implying argument_dict: :attr:`.split` = ``True``). Serves as an
        alternative to the functional group based :func:`CAT.find_substructure` function,
        which identifies the to be attached atom based on connectivity patterns
        (*i.e.* functional groups).

        In both cases the numbering of atoms starts from 1,
        following the PLAMS [1_, 2_] convention.

.. _1: https://github.com/SCM-NV/PLAMS
.. _2: https://www.scm.com/doc/plams/index.html
