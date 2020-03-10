from typing import Optional, Iterable, Any, overload, Sequence, MutableSequence, Collection, List

from rdkit import Chem
from scm.plams import Molecule

from CAT.workflow import WorkFlow
from CAT.data_handling.mol_import import read_mol
from CAT.data_handling.validate_mol import validate_mol
from CAT.attachment.ligand_opt import optimize_ligand, allign_axis
from CAT.attachment.ligand_anchoring import find_substructure


def init_multi_ligand(qd_df):
    """Initialize the multi-ligand attachment procedure."""
    workflow = WorkFlow.from_template(qd_df, name='multi_ligand')

    # Pull from the database; push unoptimized structures
    idx = None
    workflow(multi_lig, qd_df, columns=[], index=idx)


@overload
def multi_lig(qd_list: Iterable[Molecule],
              ligands: Iterable[str], dummy: Optional[Iterable[int]],
              f: None, **kwargs: Any):
    ...


@overload
def multi_lig(qd_list: Iterable[Molecule],
              ligands: Iterable[str], dummy: None,
              f: Optional[Iterable[float]], **kwargs: Any):
    ...


def multi_lig(qd_list, ligands, dummy=None, f=None, **kwargs):
    """Gogogo."""
    ligands = smiles_to_lig(list(ligands), kwargs['functional_groups'],
                            opt=kwargs['opt'], split=kwargs['split'])


def smiles_to_lig(smiles: MutableSequence[str],
                  functional_groups: Collection[Chem.Mol],
                  opt: bool = True, split: bool = True) -> List[Molecule]:
    # Convert the SMILES strings into ligands
    validate_mol(smiles, 'input_ligands')
    _ligands = read_mol(smiles)
    ligands = [find_substructure(lig, functional_groups, split)[0] for lig in _ligands]

    # Optimize the ligands
    process_mol = optimize_ligand if opt else lambda n: allign_axis(n, n.properties.dummies)
    for lig in ligands:
        process_mol(lig)
    return ligands
