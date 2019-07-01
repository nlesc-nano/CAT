"""A module designed for the calculation of Bond Dissociation Energies (BDE)."""

from itertools import chain, combinations, product
from typing import (Callable, Optional, Iterable, Tuple, Sequence, Dict, List)

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import PeriodicTable, AMSJob
from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.core.functions import (init, finish, config)
from scm.plams.core.settings import Settings

import qmflows
from data_CAT import Database

from .jobs import (job_single_point, job_geometry_opt, job_freq)
from ..utils import (get_time, type_to_string)
from ..mol_utils import (to_atnum, merge_mol)
from ..attachment.ligand_attach import rot_mol_angle
from ..properties_dataframe import PropertiesDataFrame

__all__ = ['init_bde']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
JOB_SETTINGS_BDE = ('job_settings_BDE', '')
SETTINGS1 = ('settings', 'BDE 1')
SETTINGS2 = ('settings', 'BDE 2')


def init_bde(qd_df: PropertiesDataFrame) -> None:
    """ Initialize the bond dissociation energy calculation; involves 4 distinct steps.

    * Take :math:`n` ligands (X) and another atom from the core (Y, *e.g.* Cd) and create YX*n*.
    * Given a radius :math:`r`, dissociate all possible YX*n* pairs.
    * Calculate dE: the "electronic" component of the bond dissociation energy (BDE).
    * (Optional) Calculate ddG: the thermal and entropic component of the BDE.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Unpack arguments
    path = qd_df.properties.optional.database.dirname
    overwrite = 'qd' in qd_df.properties.optional.database.overwrite
    read = 'qd' in qd_df.properties.optional.database.read
    job2 = qd_df.properties.optional.qd.dissociate.job2
    s2 = qd_df.properties.optional.qd.dissociate.s2

    # Check if the calculation has been done already
    if not overwrite and read:
        data = Database(path)
        with data.open_csv_qd(data.csv_qd, write=False) as db:
            key_ar = np.array(['BDE label', 'BDE dE', 'BDE dG', 'BDE ddG'])
            bool_ar = np.isin(key_ar, db.columns.levels[0])
            for i in db[key_ar[bool_ar]]:
                qd_df[i] = np.nan
            data.from_csv(qd_df, database='QD', get_mol=False)
        qd_df.dropna(axis='columns', how='all', inplace=True)

    # Calculate the BDEs with thermochemical corrections
    if job2 and s2:
        _bde_w_dg(qd_df)

    # Calculate the BDEs without thermochemical corrections
    else:
        _bde_wo_dg(qd_df)


def _bde_w_dg(qd_df: PropertiesDataFrame) -> None:
    """Calculate the BDEs with thermochemical corrections.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Unpack arguments
    properties = qd_df.properties
    job1 = properties.optional.qd.dissociate.job1
    job2 = properties.optional.qd.dissociate.job2
    s1 = properties.optional.qd.dissociate.s1
    s2 = properties.optional.qd.dissociate.s2
    ion = properties.optional.qd.dissociate.core_atom
    lig_count = properties.optional.qd.dissociate.lig_count
    core_index = properties.optional.qd.dissociate.core_index
    write = 'qd' in properties.optional.database.write

    # Identify previously calculated results
    try:
        has_na = qd_df[['BDE dE', 'BDE dG']].isna().all(axis='columns')
        if not has_na.any():
            return
    except KeyError:
        has_na = pd.Series(True, index=qd_df.index)

    for idx, mol in qd_df[MOL][has_na].iteritems():
        # Create XYn and all XYn-dissociated quantum dots
        xyn = get_xy2(mol, ion, lig_count)
        if not core_index:
            mol_wo_xyn = dissociate_ligand(mol, properties)
        else:
            mol_wo_xyn = dissociate_ligand2(mol, properties)

        # Construct new columns for **qd_df**
        labels = [m.properties.df_index for m in mol_wo_xyn]
        sub_idx = np.arange(len(labels)).astype(str, copy=False)
        try:
            n = qd_df['BDE label'].shape[1]
        except KeyError:
            n = 0
        if len(labels) > n:
            for i in sub_idx[n:]:
                qd_df[('BDE label', i)] = qd_df[('BDE dE', i)] = qd_df[('BDE ddG', i)] = np.nan

        # Prepare slices
        label_slice = idx, list(product(['BDE label'], sub_idx))
        dE_slice = idx, list(product(['BDE dE'], sub_idx))
        ddG_slice = idx, list(product(['BDE ddG'], sub_idx))

        # Run the BDE calculations
        init(path=mol.properties.path, folder='BDE')
        config.default_jobmanager.settings.hashing = None
        mol.properties.job_path = []
        qd_df.loc[label_slice] = labels
        qd_df.loc[dE_slice] = get_bde_dE(mol, xyn, mol_wo_xyn, job=job1, s=s1)
        qd_df.loc[ddG_slice] = get_bde_ddG(mol, xyn, mol_wo_xyn, job=job2, s=s2)
        mol.properties.job_path += xyn.properties.pop('job_path')
        for m in mol_wo_xyn:
            mol.properties.job_path += m.properties.pop('job_path')
        finish()

    qd_df['BDE dG'] = qd_df['BDE dE'] + qd_df['BDE ddG']

    job_settings = []
    for mol in qd_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    qd_df[JOB_SETTINGS_BDE] = job_settings

    # Update the database
    if write:
        with pd.option_context('mode.chained_assignment', None):
            _qd_to_db(qd_df, has_na, with_dg=True)


def _bde_wo_dg(qd_df: PropertiesDataFrame) -> None:
    """ Calculate the BDEs without thermochemical corrections.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Unpack arguments
    properties = qd_df.properties
    job1 = properties.optional.qd.dissociate.job1
    s1 = properties.optional.qd.dissociate.s1
    ion = properties.optional.qd.dissociate.core_atom
    lig_count = properties.optional.qd.dissociate.lig_count
    core_index = properties.optional.qd.dissociate.core_index
    write = 'qd' in properties.optional.database.write

    # Identify previously calculated results
    try:
        has_na = qd_df['BDE dE'].isna().all(axis='columns')
        if not has_na.any():
            return
    except KeyError:
        has_na = pd.Series(True, index=qd_df.index)

    for idx, mol in qd_df[MOL][has_na].iteritems():
        # Create XYn and all XYn-dissociated quantum dots
        xyn = get_xy2(mol, ion, lig_count)
        if not core_index:
            mol_wo_xyn = dissociate_ligand(mol)
        else:
            mol_wo_xyn = dissociate_ligand2(mol)

        # Construct new columns for **qd_df**
        labels = [m.properties.df_index for m in mol_wo_xyn]
        sub_idx = np.arange(len(labels)).astype(str, copy=False)
        try:
            n = qd_df['BDE label'].shape[1]
        except KeyError:
            n = 0
        if len(labels) > n:
            for i in sub_idx[n:]:
                qd_df[('BDE label', i)] = qd_df[('BDE dE', i)] = np.nan

        # Prepare slices
        label_slice = idx, list(product(['BDE label'], sub_idx))
        dE_slice = idx, list(product(['BDE dE'], sub_idx))

        # Run the BDE calculations
        init(path=mol.properties.path, folder='BDE')
        config.default_jobmanager.settings.hashing = None
        mol.properties.job_path = []
        qd_df.loc[label_slice] = labels
        qd_df.loc[dE_slice] = get_bde_dE(mol, xyn, mol_wo_xyn, job=job1, s=s1)
        mol.properties.job_path += xyn.properties.pop('job_path')
        for m in mol_wo_xyn:
            mol.properties.job_path += m.properties.pop('job_path')
        finish()

    job_settings = []
    for mol in qd_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    qd_df[JOB_SETTINGS_BDE] = job_settings

    # Update the database
    if write:
        with pd.option_context('mode.chained_assignment', None):
            _qd_to_db(qd_df, has_na, with_dg=False)


def _qd_to_db(qd_df: PropertiesDataFrame,
              idx: pd.Series,
              with_dg: bool = True) -> None:
    # Unpack arguments
    path = qd_df.properties.optional.database.dirname
    overwrite = 'qd' in qd_df.properties.optional.database.overwrite
    j1 = qd_df.properties.optional.qd.dissociate.job1
    s1 = qd_df.properties.optional.qd.dissociate.s1

    data = Database(path)

    qd_df.sort_index(axis='columns', inplace=True)
    kwarg = {'database': 'QD', 'overwrite': overwrite}
    if with_dg:
        j2 = qd_df.properties.optional.qd.dissociate.job2
        s2 = qd_df.properties.optional.qd.dissociate.s2
        kwarg['job_recipe'] = get_recipe(j1, s1, j2, s2)
        kwarg['columns'] = [JOB_SETTINGS_BDE, SETTINGS1, SETTINGS2]
        column_tup = ('BDE label', 'BDE dE', 'BDE ddG', 'BDE dG')
    else:
        kwarg['job_recipe'] = get_recipe(j1, s1)
        kwarg['columns'] = [JOB_SETTINGS_BDE, SETTINGS1]
        column_tup = ('BDE label', 'BDE dE')
    kwarg['columns'] += [(i, j) for i, j in qd_df.columns if i in column_tup]

    data.update_csv(qd_df[idx], **kwarg)


def get_recipe(job1: Callable,
               s1: Settings,
               job2: Optional[Callable] = None,
               s2: Optional[Callable] = None) -> Settings:
    """Return the a dictionary with job types and job settings."""
    ret = Settings()
    value1 = qmflows.singlepoint['specific'][type_to_string(job1)].copy()
    value1.update(s1)
    ret['BDE 1'] = {'key': job1, 'value': value1}

    if job2 is not None and s2 is not None:
        value2 = qmflows.freq['specific'][type_to_string(job2)].copy()
        value2.update(s2)
        ret['BDE 2'] = {'key': job2, 'value': value2}

    return ret


def get_bde_dE(tot: Molecule,
               lig: Molecule,
               core: Iterable[Molecule],
               job: Callable,
               s: Settings) -> np.ndarray:
    """Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))."""
    # Optimize XYn
    if job == AMSJob:
        s_cp = s.copy()
        s_cp.input.ams.GeometryOptimization.coordinatetype = 'Cartesian'
        lig.job_geometry_opt(job, s_cp, name='BDE_geometry_optimization')
    else:
        lig.job_geometry_opt(job, s, name='BDE_geometry_optimization')

    E_lig = lig.properties.energy.E
    if E_lig is np.nan:
        print(get_time() + 'WARNING: The BDE XYn geometry optimization failed, skipping further \
              jobs')
        return np.full(len(core), np.nan)

    # Perform a single point on the full quantum dot
    tot.job_single_point(job, s, name='BDE_single_point')
    E_tot = tot.properties.energy.E
    if E_tot is np.nan:
        print(get_time() + 'WARNING: The BDE quantum dot single point failed, \
              skipping further jobs')
        return np.full(len(core), np.nan)

    # Perform a single point on the quantum dot(s) - XYn
    for mol in core:
        mol.job_single_point(job, s, name='BDE_single_point')
    E_core = np.array([mol.properties.energy.E for mol in core])

    # Calculate and return dE
    dE = (E_lig + E_core) - E_tot
    return dE


def get_bde_ddG(tot: Molecule,
                lig: Molecule,
                core: Iterable[Molecule],
                job: Callable,
                s: Settings) -> np.ndarray:
    """Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))."""
    # Optimize XYn
    s.input.ams.Constraints.Atom = lig.properties.indices
    lig.job_freq(job, s, name='BDE_frequency_analysis')
    G_lig = lig.properties.energy.G
    E_lig = lig.properties.energy.E
    if np.nan in (E_lig, G_lig):
        print(get_time() + 'WARNING: The BDE XYn geometry optimization + freq analysis failed, \
              skipping further jobs')
        return np.full(len(core), np.nan)

    # Optimize the full quantum dot
    s.input.ams.Constraints.Atom = tot.properties.indices
    tot.job_freq(job, s, name='BDE_frequency_analysis')
    G_tot = tot.properties.energy.G
    E_tot = tot.properties.energy.E
    if np.nan in (E_tot, G_tot):
        print(get_time() + 'WARNING: The BDE quantum dot geometry optimization + freq analysis \
              failed, skipping further jobs')
        return np.full(len(core), np.nan)

    # Optimize the quantum dot(s) - XYn
    for mol in core:
        s.input.ams.Constraints.Atom = mol.properties.indices
        mol.job_freq(job, s, name='BDE_frequency_analysis')
    G_core = np.array([mol.properties.energy.G for mol in core])
    E_core = np.array([mol.properties.energy.E for mol in core])

    # Calculate and return dG and ddG
    dG = (G_lig + G_core) - G_tot
    dE = (E_lig + E_core) - E_tot
    ddG = dG - dE
    return ddG


def get_xy2(mol: Molecule,
            ion: str = 'Cd',
            lig_count: int = 2) -> Molecule:
    """Takes a quantum dot with ligands (Y) and an ion (X) and turns it into YXn.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule containing with ligands (Y).

    ion : str
        An atomic symbol (X).

    lig_count : int
        The number of ligand (*n*) in XYn.

    Returns
    -------
    plams.Molecule:
        A new XYn molecule.

    """
    def get_anchor(mol: Molecule) -> Tuple[int, Atom]:
        """Return an index and atom if marked with the properties.anchor attribute."""
        for i, at in enumerate(mol.atoms):
            if at.properties.anchor:
                return i, at

    def get_ligand(mol: Molecule) -> Molecule:
        """Extract a single ligand from **mol**."""
        at_list = []
        res = mol.atoms[-1].properties.pdb_info.ResidueNumber
        for at in reversed(mol.atoms):
            if at.properties.pdb_info.ResidueNumber == res:
                at_list.append(at)
            else:
                ret = Molecule()
                ret.atoms = at_list
                ret.bonds = list(set(chain.from_iterable(at.bonds for at in at_list)))
                return ret.copy()

    # Translate the ligands to their final position
    lig1 = get_ligand(mol)
    lig2 = lig1.copy()
    idx1, anchor1 = get_anchor(lig1)
    idx2, anchor2 = get_anchor(lig2)

    # Return a the ligand without the ion
    if ion is None:
        lig1.properties.name = 'XYn'
        lig1.properties.path = mol.properties.path
        lig1.properties.indices = [idx1]
        return lig1

    radius = anchor1.radius + PeriodicTable.get_radius(ion)
    target = np.array([radius, 0.0, 0.0])
    lig1.translate(anchor1.vector_to(target))
    lig2.translate(anchor2.vector_to(-target))

    # Define vectors for the ligand rotation
    vec1_1 = np.array(anchor1.vector_to(lig1.get_center_of_mass()))
    vec2_1 = -1 * np.array(anchor1.vector_to(np.zeros(3)))
    vec1_2 = np.array(anchor2.vector_to(lig2.get_center_of_mass()))
    vec2_2 = -1 * np.array(anchor2.vector_to(np.zeros(3)))

    # Rotate the ligands
    lig1_ar = rot_mol_angle(lig1, vec1_1, vec2_1, idx=idx1, atoms_other=anchor1, bond_length=False)
    lig2_ar = rot_mol_angle(lig2, vec1_2, vec2_2, idx=idx2, atoms_other=anchor2, bond_length=False)
    lig1.from_array(lig1_ar)
    lig2.from_array(lig2_ar)

    # Construct the CdX2 molecule
    CdX2 = Molecule()
    CdX2.add_atom(Atom(atnum=to_atnum(ion)))
    CdX2.merge_mol([lig1, lig2])
    CdX2.properties.name = 'XYn'
    CdX2.properties.path = mol.properties.path
    CdX2.properties.indices = [1, 1 + idx1, 2 + len(lig2) + idx2]
    CdX2[1].properties.charge = 0 - sum([at.properties.charge for at in CdX2.atoms[1:]])
    CdX2.properties.job_path = []

    return CdX2


def dissociate_ligand(mol: Molecule,
                      arg: Settings) -> List[Molecule]:
    """Create all XYn dissociated quantum dots.

    Parameter
    ---------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    arg: |plams.Settings|_
        A settings object containing all (optional) arguments.

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of XYn dissociated quantum dots.

    """
    # Unpack arguments
    atnum = arg.optional.qd.dissociate.core_atom
    l_count = arg.optional.qd.dissociate.lig_count
    cc_dist = arg.optional.qd.dissociate.core_core_dist
    lc_dist = arg.optional.qd.dissociate.lig_core_dist
    top_dict = arg.optional.qd.dissociate.topology

    # Convert **mol** to an XYZ array
    mol.set_atoms_id()
    xyz_array = mol.as_array()

    # Create a nested list of atoms,
    # each nested element containing all atoms with a given residue number
    res_list = []
    for at in mol:
        try:
            res_list[at.properties.pdb_info.ResidueNumber - 1].append(at)
        except IndexError:
            res_list.append([at])

    # Create a list of all core indices and ligand anchor indices
    idx_c_old = np.array([j for j, at in enumerate(res_list[0]) if at.atnum == atnum])
    idx_c, topology = filter_core(xyz_array, idx_c_old, top_dict, cc_dist)
    idx_l = np.array([i for i in mol.properties.indices if
                      mol[i].properties.pdb_info.ResidueName == 'LIG']) - 1

    # Mark the core atoms with their topologies
    for i, top in zip(idx_c_old, topology):
        mol[int(i+1)].properties.topology = top

    # Create a dictionary with core indices as keys and all combinations of 2 ligands as values
    xy = filter_lig_core(xyz_array, idx_l, idx_c, lc_dist, l_count)
    combinations_dict = get_lig_core_combinations(xy, res_list, l_count)

    # Create and return new molecules
    indices = [at.id for at in res_list[0][:-l_count]]
    indices += (idx_l[:-l_count] + 1).tolist()
    return remove_ligands(mol, combinations_dict, indices)


def dissociate_ligand2(mol: Molecule,
                       arg: Settings) -> List[Molecule]:
    """Create all XYn dissociated quantum dots.

    Parameter
    ---------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    arg: |plams.Settings|_
        A settings object containing all (optional) arguments.

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of XYn dissociated quantum dots.

    """
    # Unpack arguments
    l_count = arg.optional.qd.dissociate.lig_count
    cc_dist = arg.optional.qd.dissociate.core_core_dist
    idx_c_old = np.array(arg.optional.qd.dissociate.core_index) - 1
    top_dict = arg.optional.qd.dissociate.topology

    # Convert **mol** to an XYZ array
    mol.set_atoms_id()
    xyz_array = mol.as_array()

    # Create a nested list of atoms,
    # each nested element containing all atoms with a given residue number
    res_list = []
    for at in mol:
        try:
            res_list[at.properties.pdb_info.ResidueNumber - 1].append(at)
        except IndexError:
            res_list.append([at])

    # Create a list of all core indices and ligand anchor indices
    _, topology = filter_core(xyz_array, idx_c_old, top_dict, cc_dist)
    idx_l = np.array([i for i in mol.properties.indices if
                      mol[i].properties.pdb_info.ResidueName == 'LIG']) - 1

    # Mark the core atoms with their topologies
    for i, top in zip(idx_c_old, topology):
        mol[int(i+1)].properties.topology = top

    # Create a dictionary with core indices as keys and all combinations of 2 ligands as values
    xy = filter_lig_core2(xyz_array, idx_l, idx_c_old, l_count)
    combinations_dict = get_lig_core_combinations(xy, res_list, l_count)

    # Create and return new molecules
    indices = [at.id for at in res_list[0][:-l_count]]
    indices += (idx_l[:-l_count] + 1).tolist()
    return remove_ligands(mol, combinations_dict, indices)


def filter_lig_core2(xyz_array: np.ndarray,
                     idx_lig: Sequence[int],
                     idx_core: Sequence[int],
                     lig_count: int = 2) -> np.ndarray:
    """Create and return the indices of all possible ligand/core pairs.

    Parameters
    ----------
    xyz_array : :math:`n*3` |np.ndarray|_ [|np.float64|_]
        An array with the cartesian coordinates of a molecule with *n* atoms.

    idx_lig : |np.ndarray|_ [|np.int64|_]
        An array of all ligand anchor atoms (Y).

    idx_core : |np.ndarray|_ [|np.int64|_]
        An array of all core atoms (X).

    max_dist : float
        The maximum distance for considering XYn pairs.

    lig_count : int
        The number of ligand (*n*) in XYn.

    Returns
    -------
    :math:`m*2` |np.ndarray|_ [|np.int64|_]
        An array with the indices of all :math:`m` valid  ligand/core pairs.

    """
    dist = cdist(xyz_array[idx_lig], xyz_array[idx_core])
    xy = []
    for _ in range(lig_count):
        xy.append(np.array(np.where(dist == np.nanmin(dist, axis=0))))
        dist[xy[-1][0], xy[-1][1]] = np.nan
    xy = np.hstack(xy)
    xy = xy[[1, 0]]
    xy = xy[:, xy.argsort(axis=1)[0]]

    bincount = np.bincount(xy[0])
    xy = xy[:, [i for i, j in enumerate(xy[0]) if bincount[j] >= lig_count]]
    xy[0] = idx_core[xy[0]]
    xy[1] += 1
    return xy


def remove_ligands(mol: Molecule,
                   combinations_dict: dict,
                   indices: Sequence[int]) -> List[Molecule]:
    """ """
    ret = []
    for core in combinations_dict:
        for lig in combinations_dict[core]:
            mol_tmp = mol.copy()

            mol_tmp.properties = Settings()
            mol_tmp.properties.core_topology = str(mol[core].properties.topology) + '_' + str(core)
            mol_tmp.properties.lig_residue = sorted([mol[i[0]].properties.pdb_info.ResidueNumber
                                                     for i in lig])
            mol_tmp.properties.df_index = mol_tmp.properties.core_topology
            mol_tmp.properties.df_index += ' '.join(str(i) for i in mol_tmp.properties.lig_residue)

            delete_idx = sorted([core] + list(chain.from_iterable(lig)), reverse=True)
            for i in delete_idx:
                mol_tmp.delete_atom(mol_tmp[i])
            mol_tmp.properties.indices = indices
            mol_tmp.properties.job_path = []
            ret.append(mol_tmp)
    return ret


def filter_core(xyz_array: np.ndarray,
                idx: np.ndarray,
                topology_dict: Dict[int, str] = {6: 'vertice', 7: 'edge', 9: 'face'},
                max_dist: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """Find all atoms (**idx**) in **xyz_array** which are exposed to the surface.

    A topology is assigned to aforementioned atoms based on the number of neighbouring atoms.

    Parameters
    ----------
    xyz_array : :math:`n*3` |np.ndarray|_ [|np.float64|_]
        An array with the cartesian coordinates of a molecule with :math:`n` atoms.

    idx : |np.ndarray|_ [|np.int64|_]
        An array of atomic indices in **xyz_array**.

    topology_dict : |dict|_ [|int|_, |str|_]
        A dictionary which maps the number of neighbours (per atom) to a user-specified topology.

    max_dist : float
        The radius (Angstrom) for determining if an atom counts as a neighbour or not.

    Returns
    -------
    |np.ndarray|_ [|np.int64|_] and |np.ndarray|_ [|np.int64|_]
        The indices of all atoms in **xyz_array[idx]** exposed to the surface and
        the topology of atoms in **xyz_array[idx]**.

    """
    # Create a distance matrix and find all elements with a distance smaller than **max_dist**
    dist = cdist(xyz_array[idx], xyz_array[idx])
    np.fill_diagonal(dist, max_dist)
    xy = np.array(np.where(dist <= max_dist))
    bincount = np.bincount(xy[0], minlength=len(idx))

    # Slice xyz_array, creating arrays of reference atoms and neighbouring atoms
    x = xyz_array[idx]
    y = xyz_array[idx[xy[1]]]

    # Calculate the length of a vector from a reference atom to the mean position of its neighbours
    # A vector length close to 0.0 implies that a reference atom is surrounded by neighbours in
    # a more or less spherical pattern (i.e. the reference atom is in the bulk, not on the surface)
    vec_length = np.empty((bincount.shape[0], 3), dtype=float)
    k = 0
    for i, j in enumerate(bincount):
        vec_length[i] = x[i] - np.average(y[k:k+j], axis=0)
        k += j
    vec_length = np.linalg.norm(vec_length, axis=1)
    return idx[np.where(vec_length > 0.5)[0]], get_topology(bincount, topology_dict)


def get_topology(bincount: Iterable[int],
                 topology_dict: Dict[int, str] = {6: 'vertice', 7: 'edge', 9: 'face'}) -> List[str]:
    """Translate the number of neighbouring atoms (**bincount**) into a list of topologies.

    If a specific number of neighbours (*i*) is absent from **topology_dict** then that particular
    element is set to a generic str(*i*) + '_neighbours'.

    Parameters
    ----------
    bincount : :math:`n` |np.ndarray|_ [|np.int64|_]
        An array with the number of neighbours per atom for a total of :math:`n` atoms.

    topology_dict : |dict|_ [|int|_, |str|_]
        A dictionary which maps the number of neighbours (per atom) to a user-specified topology.

    Returns
    -------
    :math:`n` |list|_
        A list of topologies for all :math:`n` atoms in **bincount**.

    """
    if isinstance(topology_dict, Settings):
        topology_dict = topology_dict.as_dict()
    ret = []
    for i in bincount:
        try:
            ret.append(topology_dict[i])
        except KeyError:
            ret.append(str(i) + '_neighbours')
    return ret


def filter_lig_core(xyz_array: np.ndarray,
                    idx_lig: Sequence[int],
                    idx_core: Sequence[int],
                    max_dist: float = 5.0,
                    lig_count: int = 2) -> np.ndarray:
    """Create and return the indices of all possible ligand/core atom pairs.

    Ligand/core atom pair construction is limited to a given radius (**max_dist**).

    Parameters
    ----------
    xyz_array : :math:`n*3` |np.ndarray|_ [|np.float64|_]
        An array with the cartesian coordinates of a molecule with *n* atoms.

    idx_lig : |np.ndarray|_ [|np.int64|_]
        An array of all ligand anchor atoms (Y).

    idx_core : |np.ndarray|_ [|np.int64|_]
        An array of all core atoms (X).

    max_dist : float
        The maximum distance for considering XYn pairs.

    lig_count : int
        The number of ligand (*n*) in XYn.

    Returns
    -------
    :math:`m*2` |np.ndarray|_ [|np.int64|_]
        An array with the indices of all :math:`m` valid (as determined by **max_diist**)
        ligand/core pairs.

    """
    dist = cdist(xyz_array[idx_core], xyz_array[idx_lig])
    xy = np.array(np.where(dist <= max_dist))
    bincount = np.bincount(xy[0])
    xy = xy[:, [i for i, j in enumerate(xy[0]) if bincount[j] >= lig_count]]
    xy[0] = idx_core[xy[0]]
    xy[1] += 1
    return xy


def get_lig_core_combinations(xy: np.ndarray,
                              res_list: Sequence[Sequence[Atom]],
                              lig_count: int = 2) -> dict:
    """Given an array of indices (**xy**) and a nested list of atoms **res_list**.

    Parameters
    ----------
    xy : :math:`m*2` |np.ndarray|_ [|np.int64|_]
        An array with the indices of all *m* core/ligand pairs.

    res_list : |list|_ [|tuple|_ [|plams.Atom|_]]
        A list of PLAMS atoms, each nested tuple representing all atoms within a given residue.

    lig_count : int
        The number of ligand (*n*) in XYn.

    Returns
    -------
    |dict|_
        A dictionary with core/ligand combinations.

    """
    dict_ = {}
    for core, lig in xy.T:
        try:
            dict_[res_list[0][core].id].append([at.id for at in res_list[lig]])
        except KeyError:
            dict_[res_list[0][core].id] = [[at.id for at in res_list[lig]]]
    return {k: combinations(v, lig_count) for k, v in dict_.items()}
