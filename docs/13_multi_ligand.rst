.. _Multi-ligand:

Multi-ligand attachment
=======================

.. attribute:: optional.qd.multi_ligand
    :noindex:

    All settings related to the multi-ligand attachment procedure.

    Example:

    .. code:: yaml

        optional:
            qd:
                multi_ligand:
                    ligands:
                        - OCCC
                        - OCCCCCCC
                        - OCCCCCCCCCCCC
                    anchor:
                        - F
                        - Br
                        - I

|

    .. attribute:: optional.qd.multi_ligand.ligands

        :Parameter:     * **Type** - :class:`list` [:class:`str`]

        SMILES strings of to-be attached ligands.

        Note that these ligands will be attached *in addition* to whichever ligands are
        specified in :ref:`Input Cores and Ligands`.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.qd.multi_ligand.anchor

        :Parameter:     * **Type** - :class:`list` [:class:`str` or :class:`int`]

        Atomic number of symbol of the core anchor atoms.

        The first anchor atom will be assigned to the first ligand in
        :attr:`multi_ligand.ligands<optional.qd.multi_ligand.ligands>`, the second anchor atom
        to the second ligand, *etc.*.
        The list's length should consequently be of the same length as
        :attr:`multi_ligand.ligands<optional.qd.multi_ligand.ligands>`.

        Works analogous to :attr:`optional.core.anchor`.

        This optiona can alternatively be provided as ``optional.qd.multi_ligand.dummy``.

        .. note::
            This argument has no value be default and must thus be provided by the user.
