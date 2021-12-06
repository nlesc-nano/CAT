"""A module for attaching multiple non-unique ligands to a single quantum dot."""

import os
from typing import Iterable, Any, overload, Sequence, MutableSequence, List, Union

import numpy as np
import pandas as pd

from scm.plams import Molecule, MoleculeError

from .utils import AnchorTup, get_formula
from .workflows import WorkFlow
from .mol_utils import to_symbol
from .data_handling import mol_to_file
from .data_handling.mol_import import read_mol
from .data_handling.validate_mol import validate_mol
from .attachment.ligand_opt import optimize_ligand, allign_axis
from .attachment.ligand_anchoring import find_substructure
from .attachment.ligand_attach import ligand_to_qd

__all__ = ['init_multi_ligand']


def init_multi_ligand(qd_df):
    """Initialize the multi-ligand attachment procedure."""
    workflow = WorkFlow.from_template(qd_df, name='multi_ligand')

    if workflow.anchor is not None:
        sequence = [to_symbol(i) for i in workflow.anchor]
    elif workflow.f is not None:
        sequence = [str(i) for i in workflow.f]
    else:
        raise TypeError("'workflow.f' and 'workflow.anchor' cannot be both 'None'")

    columns_iter1 = ('/'.join(item for item in sequence[:i]) for i in range(1, 1+len(sequence)))
    columns_iter2 = (('multi ligand', i) for i in columns_iter1)
    columns = pd.MultiIndex.from_tuples(columns_iter2, names=qd_df.columns.names)
    for k in columns:
        qd_df[k] = None

    workflow(multi_lig, qd_df, columns=columns, no_loc=True)

    if workflow.mol_format is not None:
        path = workflow.path
        mol_format = workflow.mol_format
        mol_ar_flat = qd_df[columns].values.ravel()
        mol_to_file(mol_ar_flat, path, mol_format=mol_format)


@overload
def multi_lig(qd_series: pd.Series, ligands: Iterable[str],
              anchor: Sequence[Union[str, int]], f: None,
              **kwargs: Any) -> pd.DataFrame:
    ...
@overload  # noqa: E302
def multi_lig(qd_series: pd.Series, ligands: Iterable[str],
              anchor: None, f: Sequence[float],
              **kwargs: Any) -> pd.DataFrame:
    ...
def multi_lig(qd_series, ligands, anchor=None, f=None, **kwargs):  # noqa: E302
    """Attach multiple non-unique **ligands** to each qd in **qd_series**."""
    # Read and parse the SMILES strings
    ligands = smiles_to_lig(list(ligands),
                            functional_groups=kwargs['functional_groups'],
                            opt=kwargs['opt'],
                            split=kwargs['split'])

    mol_format = kwargs.get('mol_format')
    _path = kwargs.get('path')
    if mol_format is not _path is not None:
        path = os.path.join(str(_path).rsplit(os.sep, 1)[0], 'ligand')
        mol_to_file(ligands, path, mol_format=mol_format)

    if f is not None:
        raise NotImplementedError("'f != None' is not yet implemented")

    if anchor is not None:
        return _multi_lig_anchor(qd_series, ligands, kwargs['path'], anchor, kwargs['allignment'])
    elif f is not None:
        return [[NotImplemented]]
    else:
        raise TypeError("'f' and 'anchor' cannot be both 'None'")


def _multi_lig_anchor(qd_series, ligands, path, anchor, allignment) -> np.ndarray:
    """Gogogo."""
    ret = np.empty((len(ligands), len(qd_series)), dtype=object)
    for i, qd in enumerate(qd_series):
        qd = qd.copy()

        for j, (ligand, atnum) in enumerate(zip(ligands, anchor)):
            qd.set_atoms_id()
            try:
                atoms = [at for at in qd if at.atnum == atnum]
                assert atoms
            except AssertionError as ex:
                raise MoleculeError(f'Failed to identify {to_symbol(atnum)!r} in '
                                    f'{get_formula(qd)!r}') from ex

            coords = Molecule.as_array(None, atom_subset=atoms)
            qd.properties.dummies = np.array(coords, ndmin=2, dtype=float)
            qd = ligand_to_qd(qd, ligand, path=path,
                              allignment=allignment,
                              idx_subset=qd.properties.indices)
            ret[j, i] = qd
            for at in reversed(atoms):
                qd.delete_atom(qd[at.id])
            qd.unset_atoms_id()
    return ret


def _multi_lig_f(qd_series, ligands, path, f, **kwargs):
    raise NotImplementedError


def smiles_to_lig(smiles: MutableSequence[str],
                  functional_groups: Iterable[AnchorTup],
                  opt: bool = True, split: bool = True) -> List[Molecule]:
    """Parse and convert all **smiles** strings into Molecules."""
    # Convert the SMILES strings into ligands
    validate_mol(smiles, 'input_ligands')
    ligands = [find_substructure(lig, functional_groups, split)[0] for lig in read_mol(smiles)]

    # Optimize the ligands
    process_mol = optimize_ligand if opt else lambda n: allign_axis(n)
    for lig in ligands:
        process_mol(lig)
    return ligands
