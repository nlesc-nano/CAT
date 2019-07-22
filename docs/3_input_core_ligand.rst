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

    .. warning:
        This option works reliably only for geometries representing complete molecules.
        If some atoms are missing (for example, a protein without hydrogens) the resulting set of bonds
        would usually contain more bonds or bonds with higher order than expected.


.. attribute:: .column

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The column containing the to be imported molecules.
    Relevant when importing structures from .txt and .xlsx files with
    multiple columns.
    Relevant for .txt and .csv files.
    Numbering starts from 0.


.. attribute:: .row

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The first row in a column which contains a molecule.
    Useful for when, for example, the very first row contains the title of
    aforementioned row, in which case row = 1 would be a sensible choice.
    Relevant for .txt and .csv files.
    Numbering starts from 0.

.. attribute:: .indices

    :Parameter:     * **Type** - :class:`int` or :class:`tuple` [:class:`int`]
                    * **Default value** – ``None``

    The behaviour of this argument depends on whether it is passed to a molecule
    in :attr:`input_cores` or :attr:`input_ligands`:

    .. attribute:: input_cores

        Manually specify the atomic index of one ore more atom(s) in the core that
        will be replaced with ligands. If left empty, all atoms of a user-specified
        element (see :attr:`optional.cores.dummy`) will be replaced with
        ligands.

    .. attribute:: input_ligands

        Manually specify the atomic index of the ligand atom that will be attached
        to core (implying argument_dict: :attr:`optional.ligand.split` = ``False``).
        If two atomic indices are provided (*e.g.* ``(1, 2)``), the bond between atoms ``1`` and
        [``2``] will be broken and the remaining molecule containing atom ``2`` is attached to the core,
        (implying argument_dict: :attr:`.split` = ``True``).
        Serves as an alternative to the functional group based :func:`CAT.find_substructure` function,
        which identifies the to be attached atom based on connectivity patterns
        (*i.e.* functional groups).

    .. note::
        Atom numbering follows the PLAMS [1_, 2_] convention of starting from 1 rather than 0.

.. _1: https://github.com/SCM-NV/PLAMS
.. _2: https://www.scm.com/doc/plams/index.html
