""" A module designed for the calculation of Bond Dissociation Energies (BDE). """

__all__ = ['init_bde']

from itertools import chain, combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.core.functions import (init, finish, config)
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob

from .jobs import (job_single_point, job_geometry_opt, job_freq)
from .. import utils as CAT
from ..utils import get_time
from ..mol_utils import (to_atnum, merge_mol)
from ..attachment.ligand_attach import rot_mol_angle
from ..data_handling.database import (property_to_database, get_empty_columns, _anchor_to_idx)


def init_bde(mol_list, arg):
    """ Initialize the bond dissociation energy calculation; involves 4 distinct steps:
    1.  Take two ligands X and another atom from the core Y (e.g. Cd) and create YX2.
    2.  Create all n*2*(n-1) possible molecules where YX2 is dissociated.
    3.  Calculate dE: the "electronic" component of the bond dissociation energy (BDE).
    4.  Calculate ddG: the thermal and entropic component of the BDE.

    mol_list <list> [<plams.Molecule>]: A list of PLAMS molecule.
    job1 <type> & s1 <Settings>: A type object of a job and its settings; used in step 3.
    job2 <type> & s2 <Settings>: A type object of a job and its settings; used in step 4.
    return <pd.DataFrame>: A pandas dataframe with ligand residue numbers, Cd topology and BDEs.
    """
    # Prepare the job settings
    job_recipe = arg.optional.qd.dissociate

    # Check if the calculation has been done already
    if 'qd' not in arg.optional.database.overwrite:
        empty_columns = get_empty_columns('BDE dE', arg, database='qd')
        mol_list = [mol for mol in mol_list if
                    (mol.properties.core, mol.properties.core_anchor,
                     mol.properties.ligand, mol.properties.ligand_anchor) in empty_columns]
        if not mol_list:
            return

    # Prepare the dataframe
    df = _get_bde_df()
    for mol in mol_list:
        # Ready YX2 and the YX2 dissociated quantum dots
        lig = get_cdx2(mol)
        core = dissociate_ligand(mol, arg)

        # Construct a Series
        index = [i.properties.df_index for i in core]
        series = _get_bde_series(index, mol, job_recipe.job2)

        # Fill the series with energies
        init(path=mol.properties.path, folder='BDE')
        config.default_jobmanager.settings.hashing = None
        series['BDE dE'] = get_bde_dE(mol, lig, core, job=job_recipe.job1, s=job_recipe.s1)
        if job_recipe.job2 and job_recipe.s2:
            series['BDE ddG'] = get_bde_ddG(mol, lig, core, job=job_recipe.job2, s=job_recipe.s2)
            series['BDE dG'] = series['BDE dE'] + series['BDE ddG']
        finish()

        # Update the indices of df
        for i in series.index:
            if i not in df.index:
                df.loc[i, :] = None

        # Update the values of df
        df[series.name] = series

    # Export the BDE results to the database
    del df[(None, None, None, None)]
    if 'qd' in arg.optional.database.write:
        property_to_database(df, arg, database='qd', prop='bde')


def _get_bde_df():
    """ Return an empty dataframe for init_bde(). """
    idx_names = ['index', 'sub index']
    idx = pd.MultiIndex(levels=[[], []], codes=[[], []], names=idx_names)
    column_names = ['core', 'core_anchor', 'ligand smiles', 'ligand anchor']
    columns = pd.MultiIndex.from_tuples([(None, None, None, None)], names=column_names)
    return pd.DataFrame(index=idx, columns=columns)


def _get_bde_series(index, mol, job2=True):
    """ Return an empty series for the for loop in init_bde(). """
    name = (mol.properties.core, mol.properties.core_anchor,
            mol.properties.ligand, mol.properties.ligand_anchor)
    if job2:
        super_index = ['BDE dE', 'BDE ddG', 'BDE dG']
    else:
        super_index = ['BDE dE']
    index = [(i, j) for j in index for i in super_index]
    index.sort()
    index.append(('BDE settings1', ''))
    if job2:
        index.append(('BDE settings2', ''))
    idx = pd.MultiIndex.from_tuples(index, names=['index', 'sub index'])
    return pd.Series(None, index=idx, name=name)


def get_bde_dE(tot, lig, core, job=None, s=None):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
    # Switch to default settings if no job & s are <None>
    if job is None and s is None:
        job = AMSJob
        s = CAT.get_template('qd.yaml')['MOPAC']
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

    # Optimize XYn
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


def get_bde_ddG(tot, lig, core, job=None, s=None):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
    # Switch to default settings if no job & s are <None>
    if job is None and s is None:
        job = AMSJob
        s = CAT.get_template('qd.yaml')['UFF']
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

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


def get_cdx2(mol, ion='Cd'):
    """ Takes a quantum dot with ligands (X) and an ion (Y) and turns it into YX2.
    Returns the total energy of YX2 at the MOPAC level of theory. """
    def get_anchor(mol):
        """ Return an index and atom if marked with the properties.anchor attribute """
        for i, at in enumerate(mol.atoms):
            if at.properties.anchor:
                return i, at

    def get_ligand(mol):
        """ Extract a single ligand from *mol*. """
        at_list = []
        res = mol.atoms[-1].properties.pdb_info.ResidueNumber
        for at in reversed(mol.atoms):
            if at.properties.pdb_info.ResidueNumber == res:
                at_list.append(at)
            else:
                ret = Molecule()
                ret.atoms = at_list
                ret.bonds = set(chain.from_iterable(at.bonds for at in at_list))
                return ret.copy()

    # Translate the ligands to their final position
    lig1 = get_ligand(mol)
    lig2 = lig1.copy()
    idx1, anchor1 = get_anchor(lig1)
    idx2, anchor2 = get_anchor(lig2)
    target = np.array([2.2, 0.0, 0.0])
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
    CdX2.properties.name = 'YX2'
    CdX2.properties.path = mol.properties.path
    CdX2.properties.indices = [1, 1 + idx1, 2 + len(lig2) + idx2]

    return CdX2


def dissociate_ligand(mol, arg):
    """ """
    # Unpack arguments
    atnum = arg.optional.qd.dissociate.core_atom
    lig_count = arg.optional.qd.dissociate.lig_count
    core_core_dist = arg.optional.qd.dissociate.core_core_dist
    lig_core_dist = arg.optional.qd.dissociate.lig_core_dist
    topology_dict = arg.optional.qd.dissociate.topology

    # Convert mol
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
    idx_core_old = np.array([j for j, at in enumerate(res_list[0]) if at.atnum == atnum])
    idx_core, topology = filter_core(xyz_array, idx_core_old, topology_dict, core_core_dist)
    i, j = len(res_list[0]), len(res_list[1])
    k = _anchor_to_idx(mol.properties.ligand_anchor)
    idx_lig = np.arange(i + k, i + k + j * len(res_list[1:]), j) - 1

    # Mark the core atoms with their topologies
    for i, top in zip(idx_core_old, topology):
        mol[int(i+1)].properties.topology = top

    # Create a dictionary with core indices as keys and all combinations of 2 ligands as values
    xy = filter_lig_core(xyz_array, idx_lig, idx_core, lig_core_dist)
    combinations_dict = get_lig_core_combinations(xy, res_list, lig_count)

    # Create and return new molecules
    indices = [at.id for at in res_list[0][:-lig_count]] + (idx_lig[:-lig_count] + 1).tolist()
    return remove_ligands(mol, combinations_dict, indices, idx_lig)


def remove_ligands(mol, combinations_dict, indices, idx_lig):
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
            mol_tmp.properties.df_index += ''.join([' ' + str(i)
                                                    for i in mol_tmp.properties.lig_residue])
            delete_idx = sorted([core + 1] + list(chain.from_iterable(lig)), reverse=True)
            for i in delete_idx:
                mol_tmp.delete_atom(mol_tmp[i])
            mol_tmp.properties.indices = indices
            ret.append(mol_tmp)
    return ret


def filter_core(xyz_array, idx, topology_dict, max_dist=5.0, ret_threshold=0.5):
    """ """
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
    # a more or less spherical pattern (i.e. it is not exposed to the surface)
    vec_length = np.empty((bincount.shape[0], 3), dtype=float)
    k = 0
    for i, j in enumerate(bincount):
        vec_length[i] = x[i] - np.average(y[k:k+j], axis=0)
        k += j
    vec_length = np.linalg.norm(vec_length, axis=1)

    return idx[np.where(vec_length > ret_threshold)[0]], get_topology(bincount, topology_dict)


def get_topology(bincount, topology_dict={6: 'vertice', 7: 'edge', 9: 'face'}):
    """ """
    if isinstance(topology_dict, Settings):
        topology_dict = topology_dict.as_dict()
    ret = []
    for i in bincount:
        try:
            ret.append(topology_dict[i])
        except KeyError:
            ret.append(str(i) + '_neighbours')
    return ret


def filter_lig_core(xyz_array, idx_lig, idx_core, max_dist=5.0):
    """ """
    dist = cdist(xyz_array[idx_lig], xyz_array[idx_core]).T
    xy = np.array(np.where(dist <= max_dist))
    bincount = np.bincount(xy[0])
    xy = xy[:, [i for i, j in enumerate(xy[0]) if bincount[j] >= 2]]
    xy[0] = idx_core[xy[0]]
    xy[1] += 1

    return xy


def get_lig_core_combinations(xy, res_list, lig_count=2):
    """ """
    ret = {}
    for x, y in xy.T:
        try:
            ret[res_list[0][x].id].append([at.id for at in res_list[y]])
        except KeyError:
            ret[res_list[0][x].id] = [[at.id for at in res_list[y]]]
    for i in ret:
        ret[i] = combinations(ret[i], lig_count)

    return ret
