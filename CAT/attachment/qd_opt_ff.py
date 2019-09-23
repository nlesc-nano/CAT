"""
CAT.attachment.qd_opt_ff
========================

A module designed for optimizing the combined ligand & core using with forcefields in CP2K.

Index
-----
.. currentmodule:: CAT.attachment.qd_opt_ff
.. autosummary::
    qd_opt_ff
    get_psf

API
---
.. autofunction:: qd_opt_ff
.. autofunction:: get_psf

"""

import os
from collections import abc
from typing import Container, Iterable, Union

import numpy as np
import pandas as pd

from scm.plams import Molecule, Settings

try:
    from nanoCAT.ff.cp2k_utils import set_cp2k_element
    from nanoCAT.ff.psf import PSF
    NANO_CAT: bool = True
except ImportError:
    PSF = 'PSF'
    NANO_CAT: bool = False

__all__ = ['qd_opt_ff']


def qd_opt_ff(mol: Molecule, job_recipe: Settings, name: str = 'QD_opt') -> None:
    """Alternative implementation of :func:`.qd_opt` using CP2Ks' classical forcefields.

    Performs an inplace update of **mol**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        The to-be optimized molecule.

    job_recipe : |plams.Settings|_
        A Settings instance containing all jon settings.
        Expects 4 keys: ``"job1"``, ``"job2"``, ``"s1"``, ``"s2"``.

    name : str
        The name of the job.

    See also
    --------
    :func:`CAT.attachment.qd_opt.qd_opt`
        Default workflow for optimizing molecules.

    """
    psf_name = os.path.join(mol.properties.path, mol.properties.name + '.psf')

    # Prepare the job settings
    job, s = job_recipe.job1, job_recipe.s1.copy()

    s.runscript.pre = (f'ln "{psf_name}" ./"{name}.psf"\n'
                       f'ln "{mol.properties.prm}" ./"{name}.prm"')
    s.input.force_eval.subsys.topology.conn_file_name = f'{name}.psf'
    s.input.force_eval.mm.forcefield.parm_file_name = f'{name}.prm'
    set_cp2k_element(s, mol)

    if not os.path.isfile(psf_name):
        psf = get_psf(mol, s.input.force_eval.mm.forcefield.charge)
        psf.write(psf_name)

    # Run the first job and fix broken angles
    mol.job_geometry_opt(job, s, name=name, read_template=False)
    mol.round_coords()


def get_psf(mol: Molecule, charges: Union[Settings, Iterable[Settings]]) -> PSF:
    """Construct and return a :class:`PSF` instance.

    .. _CHARGE: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/CHARGE.html

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule which will be used for constructing the :class:`PSF` instance.

    s : |list|_ [|plams.Settings|_] or |plams.Settings|_
        A list of settings constructed from the CP2K /FORCE_EVAL/MM/FORCEFIELD/CHARGE_ block.
        The settings are expected to contain the ``"charge"`` and ``"atom"`` keys.

    """
    # Construct a PSF instance
    psf = PSF()
    psf.generate_bonds(mol)
    psf.generate_angles(mol)
    psf.generate_dihedrals(mol)
    psf.generate_impropers(mol)
    psf.generate_atoms(mol)

    # Update charges based on charges which have been explictly specified by the user
    initial_charge = psf.charge.sum()
    if isinstance(charges, Settings):
        charge_dict = {charges.atom: charges.charge}
    elif isinstance(charges, abc.Iterable):
        charge_dict = {i.atom: i.charge for i in charges}
    else:
        raise TypeError(f"The parameter 'charges' is of invalid type: {repr(type(charges))}")

    for at, charge in charge_dict.items():
        psf.update_atom_charge(at, float(charge))

    # Update atomic charges in order to reset the molecular charge to its initial value
    constrain_charge(psf, initial_charge, charge_dict)
    return psf


def constrain_charge(psf: PSF, initial_charge: float, atom_set: Container[str]) -> None:
    """Set to total molecular charge of **psf** to **initial_charge**.

    Atoms in **psf** whose atomic symbol intersects with **charge_set** will *not*
    be altered.

    Parameters
    ----------
    psf : |nanoCAT.PSF|_
        A :class:`.PSF` instance with newly updated charges and/or atom types.

    initial_charge : float
        The initial charge of the system before the updating of atomic charges.

    charge_set : |Container|_ [|str|_]
        A container with atomic symbols.
        Any atom in **psf** whose atomic symbol intersects with **charge_set** will *not*
        be altered.

    """
    # Check if the molecular charge has remained unchanged
    new_charge = psf.charge.sum()
    if (abs(new_charge) - abs(initial_charge)) < 10**-8:  # i.e. 0.0 +- 10**-8
        return psf

    # Update atomic charges in order to reset the molecular charge to its initial value
    atom_subset = np.array([at not in atom_set for at in psf.atom_type])
    charge_correction = initial_charge - psf.charge[atom_subset].sum()
    charge_correction /= np.count_nonzero(atom_subset)
    with pd.option_context('mode.chained_assignment', None):
        psf.charge[atom_subset] += charge_correction
