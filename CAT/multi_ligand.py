from typing import (Iterable, Any, overload, Sequence, MutableSequence,
                    Collection, List, Generator, Union)

import numpy as np
import pandas as pd

from rdkit import Chem
from scm.plams import Molecule

from CAT.workflows import WorkFlow, MOL
from CAT.utils import group_by_values
from CAT.mol_utils import to_symbol
from CAT.data_handling import mol_to_file
from CAT.data_handling.mol_import import read_mol
from CAT.data_handling.validate_mol import validate_mol
from CAT.attachment.ligand_opt import optimize_ligand, allign_axis
from CAT.attachment.ligand_anchoring import find_substructure
from CAT.attachment.ligand_attach import ligand_to_qd


def init_multi_ligand(qd_df):
    """Initialize the multi-ligand attachment procedure."""
    workflow = WorkFlow.from_template(qd_df, name='multi_ligand')

    if workflow.dummy is not None:
        sequence = [to_symbol(i) for i in workflow.dummy]
    elif workflow.f is not None:
        sequence = [str(i) for i in workflow.f]
    else:
        raise TypeError("'workflow.f' and 'workflow.dummy' cannot be both 'None'")

    columns_iter1 = ('/'.join(item for item in sequence[:i]) for i in range(1, 1+len(sequence)))
    columns_iter2 = (('multi ligand', i) for i in columns_iter1)
    columns = pd.MultiIndex.from_tuples(columns_iter2, names=qd_df.columns.names)

    workflow(multi_lig, qd_df[MOL], columns=columns)

    if workflow.mol_format is not None:
        path = workflow.path
        mol_format = workflow.mol_format
        mol_to_file(qd_df[columns].values.ravel(), path, mol_format=mol_format)


@overload
def multi_lig(qd_series: pd.Series, ligands: Iterable[str], columns: pd.MultiIndex,
              dummy: Sequence[Union[str, int]], f: None,
              **kwargs: Any) -> pd.DataFrame:
    ...


@overload
def multi_lig(qd_series: pd.Series, ligands: Iterable[str], columns: pd.MultiIndex,
              dummy: None, f: Sequence[float],
              **kwargs: Any) -> pd.DataFrame:
    ...


def multi_lig(qd_series, ligands, columns, dummy=None, f=None, **kwargs):
    ligands = smiles_to_lig(list(ligands),
                            functional_groups=kwargs['functional_groups'],
                            opt=kwargs['opt'],
                            split=kwargs['split'])

    if f is not None:
        raise NotImplementedError("'f != None' is not yet implemented")

    if dummy is not None:
        data = _multi_lig_dummy(qd_series, ligands, kwargs['path'], dummy, kwargs['allignment'])
    elif f is not None:
        data = NotImplemented
    else:
        raise TypeError("'f' and 'dummy' cannot be both 'None'")

    return pd.DataFrame(data, index=qd_series.index, columns=columns, dtype=object)


def _multi_lig_dummy(qd_series, ligands, path, dummy, allignment
                     ) -> Generator[List[Molecule], None, None]:
    """Gogogo."""
    for qd in qd_series:
        atnum_dict = group_by_values(enumerate(at.atnum for at in qd))
        ret = []

        for ligand, atnum in zip(ligands, dummy):
            idx = atnum_dict[atnum]
            qd.properties.dummies = np.fromiter(idx, dtype=int, count=len(idx))

            qd = ligand_to_qd(qd, ligand, path=path,
                              allignment=allignment,
                              idx_subset=qd.properties.indices)
            del qd.properties.dummies
            ret.append(qd)
        yield ret


def _multi_lig_f(qd_series, ligands, path, f, **kwargs):
    raise NotImplementedError


def smiles_to_lig(smiles: MutableSequence[str],
                  functional_groups: Collection[Chem.Mol],
                  opt: bool = True, split: bool = True) -> List[Molecule]:
    """Parse and convert all **smiles** strings into Molecules."""
    # Convert the SMILES strings into ligands
    validate_mol(smiles, 'input_ligands')
    ligands = [find_substructure(lig, functional_groups, split)[0] for lig in read_mol(smiles)]

    # Optimize the ligands
    process_mol = optimize_ligand if opt else lambda n: allign_axis(n, n.properties.dummies)
    for lig in ligands:
        process_mol(lig)
    return ligands
