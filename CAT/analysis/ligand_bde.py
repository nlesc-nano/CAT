""" A module designed for the calculation of Bond Dissociation Energies (BDE). """

__all__ = ['init_bde']

from itertools import chain, combinations, product

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import PeriodicTable
from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.core.functions import (init, finish, config)
from scm.plams.core.settings import Settings

import qmflows

from .jobs import (job_single_point, job_geometry_opt, job_freq)
from ..utils import (get_time, type_to_string)
from ..mol_utils import (to_atnum, merge_mol)
from ..attachment.ligand_attach import rot_mol_angle
from ..data_handling.CAT_database import Database


def init_bde(qd_df, arg):
    """ Initialize the bond dissociation energy calculation; involves 4 distinct steps:
    1.  Take *n* ligands (X) and another atom from the core (Y, *e.g.* Cd) and create YX*n*.
    2.  Given a radius *r*, dissociate all possible YX*n* pairs.
    3.  Calculate dE: the "electronic" component of the bond dissociation energy (BDE).
    4.  (Optional) Calculate ddG: the thermal and entropic component of the BDE.

    :parameter qd_df: A dataframe of quantum dots.
    :type qd_df: |pd.DataFrame|_ (columns: |str|_, index=|str|_, values=|plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    """
    data = Database(arg.optional.database.dirname)
    overwrite = 'qd' in arg.optional.database.overwrite

    # Check if the calculation has been done already
    if not overwrite and 'qd' in arg.optional.database.read:
        with data.open_csv_qd(data.csv_qd, write=False) as db:
            try:
                for i in db[['BDE label', 'BDE dE', 'BDE dG', 'BDE ddG']]:
                    qd_df[i] = np.nan
            except KeyError:
                pass
            data.from_csv(qd_df, database='QD', get_mol=False)
        qd_df.dropna(axis='columns', how='all', inplace=True)

    # Calculate the BDEs with thermochemical corrections
    if arg.optional.qd.dissociate.job2 and arg.optional.qd.dissociate.s2:
        _bde_w_dg(qd_df, arg)

    # Calculate the BDEs without thermochemical corrections
    else:
        _bde_wo_dg(qd_df, arg)


def _bde_w_dg(qd_df, arg):
    """ Calculate the BDEs with thermochemical corrections.

    :parameter qd_df: A dataframe of quantum dots.
    :type qd_df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    """
    data = Database(arg.optional.database.dirname)
    overwrite = 'qd' in arg.optional.database.overwrite

    # Prepare the job settings
    j1, j2 = arg.optional.qd.dissociate.job1, arg.optional.qd.dissociate.job2
    s1, s2 = arg.optional.qd.dissociate.s1, arg.optional.qd.dissociate.s2

    try:
        has_na = qd_df[['BDE dE', 'BDE dG']].isna().all(axis='columns')
        if not has_na.any():
            return
    except KeyError:
        has_na = qd_df.index

    for i, mol in qd_df['mol'][has_na].iteritems():
        # Create XYn and all XYn-dissociated quantum dots
        xyn = get_xy2(mol)
        if not arg.optional.qd.dissociate.core_index:
            mol_wo_xyn = dissociate_ligand(mol, arg)
        else:
            mol_wo_xyn = dissociate_ligand2(mol, arg)

        # Construct new columns for **qd_df**
        values = [i.properties.df_index for i in mol_wo_xyn]
        n = len(values)
        sub_idx = np.arange(n)
        super_idx = ['BDE label', 'BDE dE', 'BDE ddG', 'BDE dG']
        idx = list(product(super_idx, sub_idx))
        for j in idx:
            if j not in qd_df.index:
                qd_df[j] = np.nan

        # Prepare slices
        label_slice = i, idx[0:n]
        dE_slice = i, idx[n:2*n]
        ddG_slice = i, idx[2*n:3*n]

        # Run the BDE calculations
        init(path=mol.properties.path, folder='BDE')
        config.default_jobmanager.settings.hashing = None
        qd_df.loc[label_slice] = values
        qd_df.loc[dE_slice] = get_bde_dE(mol, xyn, mol_wo_xyn, job=j1, s=s1)
        qd_df.loc[ddG_slice] = get_bde_ddG(mol, xyn, mol_wo_xyn, job=j2, s=s2)
        finish()
    qd_df['BDE dG'] = qd_df['BDE dE'] + qd_df['BDE ddG']

    # Update the database
    if 'qd' in arg.optional.database.write:
        value1 = qmflows.singlepoint['specific'][type_to_string(j1)].copy()
        value1.update(s1)
        value2 = qmflows.freq['specific'][type_to_string(j2)].copy()
        value2.update(s2)
        recipe = Settings()
        recipe['BDE 1'] = {'key': j1, 'value': value1}
        recipe['BDE 2'] = {'key': j2, 'value': value2}
        data.update_csv(qd_df, database='QD', job_recipe=recipe, overwrite=overwrite,
                        columns=[('settings', 'BDE 1'), ('settings', 'BDE 2')]+idx)


def _bde_wo_dg(qd_df, arg):
    """ Calculate the BDEs without thermochemical corrections.

    :parameter qd_df: A dataframe of quantum dots.
    :type qd_df: |pd.DataFrame|_ (columns: |str|_, index=|str|_, values=|plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    """
    data = Database(arg.optional.database.dirname)
    overwrite = 'qd' in arg.optional.database.overwrite

    # Prepare the job settings
    j1 = arg.optional.qd.dissociate.job1
    s1 = arg.optional.qd.dissociate.s1

    try:
        has_na = qd_df['BDE dE'].isna().all(axis='columns')
        if not has_na.any():
            return
    except KeyError:
        has_na = qd_df.index

    for i, mol in qd_df['mol'][has_na].iteritems():
        # Create XYn and all XYn-dissociated quantum dots
        xyn = get_xy2(mol)
        if not arg.optional.qd.dissociate.core_index:
            mol_wo_xyn = dissociate_ligand(mol, arg)
        else:
            mol_wo_xyn = dissociate_ligand2(mol, arg)

        # Construct new columns for **qd_df**
        values = [i.properties.df_index for i in mol_wo_xyn]
        n = len(values)
        sub_idx = np.arange(n)
        super_idx = ['BDE label', 'BDE dE']
        idx = list(product(super_idx, sub_idx))
        for j in idx:
            if j not in qd_df.index:
                qd_df[j] = np.nan

        # Prepare slices
        label_slice = i, idx[0:n]
        dE_slice = i, idx[n:2*n]

        # Run the BDE calculations
        init(path=mol.properties.path, folder='BDE')
        config.default_jobmanager.settings.hashing = None
        qd_df.loc[label_slice] = values
        qd_df.loc[dE_slice] = get_bde_dE(mol, xyn, mol_wo_xyn, job=j1, s=s1)
        finish()

    # Update the database
    if 'qd' in arg.optional.database.write:
        recipe = Settings()
        value = qmflows.singlepoint[type_to_string(j1)]
        value.update(s1)
        recipe['BDE 1'] = {'key': j1, 'value': value}
        data.update_csv(qd_df, database='QD', job_recipe=recipe, overwrite=overwrite,
                        columns=[('settings', 'BDE 1')]+idx)


def get_bde_dE(tot, lig, core, job=None, s=None):
    """ Calculate the bond dissociation energy: dE = dE(mopac) + (dG(uff) - dE(uff))
    """
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


def get_xy2(mol, ion='Cd', lig_count=2):
    """ Takes a quantum dot with ligands (Y) and an ion (X) and turns it into YXn.
    Returns a XYn molecule.

    :parameter mol: A PLAMS molecule containing with ligands (Y).
    :type mol: |plams.Molecule|_
    :parameter str ion: An atomic symbol (X).
    :parameter int lig_count: The number of ligand (*n*) in XYn.
    :return: A new XYn molecule.
    :rtype: plams.Molecule.
    """
    def get_anchor(mol):
        """ Return an index and atom if marked with the properties.anchor attribute """
        for i, at in enumerate(mol.atoms):
            if at.properties.anchor:
                return i, at

    def get_ligand(mol):
        """ Extract a single ligand from **mol**. """
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

    return CdX2


def dissociate_ligand(mol, arg):
    """ Create all XYn dissociated quantum dots.

    :parameter mol: A PLAMS molecule.
    :type mol: |plams.Molecule|_
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
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
    return remove_ligands(mol, combinations_dict, indices, idx_l)


def dissociate_ligand2(mol, arg):
    """ Create all XYn dissociated quantum dots.

    :parameter mol: A PLAMS molecule.
    :type mol: |plams.Molecule|_
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
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
    return remove_ligands(mol, combinations_dict, indices, idx_l)


def filter_lig_core2(xyz_array, idx_lig, idx_core, lig_count=2):
    """ Create and return the indices of all possible ligand/core pairs..

    :parameter xyz_array: An array with the cartesian coordinates of a molecule with *n* atoms.
    :type xyz_array: *n*3* |np.ndarray|_ [|np.float64|_]
    :parameter idx: An array of all ligand anchor atoms (Y).
    :type idx: |np.ndarray|_ [|np.int64|_]
    :parameter idx: An array of all core atoms (X).
    :type idx: |np.ndarray|_ [|np.int64|_]
    :parameter int lig_count: The number of ligand (*n*) in XYn.
    :return: An array with the indices of all *m* valid (as determined by **max_diist**)
        ligand/core pairs.
    :rtype: *m*2* |np.ndarray|_ [|np.int64|_].
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
            delete_idx = sorted([core] + list(chain.from_iterable(lig)), reverse=True)
            for i in delete_idx:
                mol_tmp.delete_atom(mol_tmp[i])
            mol_tmp.properties.indices = indices
            ret.append(mol_tmp)
    return ret


def filter_core(xyz_array, idx, topology_dict={6: 'vertice', 7: 'edge', 9: 'face'}, max_dist=5.0):
    """ Find all atoms (**idx**) in **xyz_array** which are exposed to the surface
    and assign a topology to aforementioned atoms based on the number of neighbouring atoms.

    :parameter xyz_array: An array with the cartesian coordinates of a molecule with *n* atoms.
    :type xyz_array: *n*3* |np.ndarray|_ [|np.float64|_]
    :parameter idx: An array of atomic indices in **xyz_array**.
    :type idx: |np.ndarray|_ [|np.int64|_]
    :parameter topology_dict: A dictionary which maps the number of neighbours (per atom) to a
        user-specified topology.
    :type topology_dict: |dict|_ (keys: |int|_)
    :parameter float max_dist: The radius (Angstrom) for determining if an atom counts as a
        neighbour or not.
    :return: The indices of all atoms in **xyz_array[idx]** exposed to the surface and
        the topology of atoms in **xyz_array[idx]**.
    :rtype: |np.ndarray|_ [|np.int64|_] and |np.ndarray|_ [|np.int64|_]
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


def get_topology(bincount, topology_dict={6: 'vertice', 7: 'edge', 9: 'face'}):
    """ Translate the number of neighbouring atoms (**bincount**) into a list of topologies.
    If a specific number of neighbours (*i*) is absent from **topology_dict** then that particular
    element is set to a generic str(*i*) + '_neighbours'.

    :parameter bincount: An array with the number of neighbours per atom for a total of *n* atoms.
    :type bincount: *n* |np.ndarray|_ [|np.int64|_]
    :parameter topology_dict: A dictionary which maps the number of neighbours (per atom) to a
        user-specified topology.
    :type topology_dict: |dict|_ (keys: |int|_)
    :return: A list of topologies for all *n* atoms in **bincount**.
    :rtype: *n* |list|_.
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


def filter_lig_core(xyz_array, idx_lig, idx_core, max_dist=5.0, lig_count=2):
    """ Create and return the indices of all possible ligand/pairs that can be constructed within a
    given radius (**max_dist**).

    :parameter xyz_array: An array with the cartesian coordinates of a molecule with *n* atoms.
    :type xyz_array: *n*3* |np.ndarray|_ [|np.float64|_]
    :parameter idx: An array of all ligand anchor atoms (Y).
    :type idx: |np.ndarray|_ [|np.int64|_]
    :parameter idx: An array of all core atoms (X).
    :type idx: |np.ndarray|_ [|np.int64|_]
    :parameter float max_dist: The maximum distance for considering XYn pairs.
    :parameter int lig_count: The number of ligand (*n*) in XYn.
    :return: An array with the indices of all *m* valid (as determined by **max_diist**)
        ligand/core pairs.
    :rtype: *m*2* |np.ndarray|_ [|np.int64|_].
    """
    dist = cdist(xyz_array[idx_core], xyz_array[idx_lig])
    xy = np.array(np.where(dist <= max_dist))
    bincount = np.bincount(xy[0])
    xy = xy[:, [i for i, j in enumerate(xy[0]) if bincount[j] >= lig_count]]
    xy[0] = idx_core[xy[0]]
    xy[1] += 1
    return xy


def get_lig_core_combinations(xy, res_list, lig_count=2):
    """ Given an array of indices (**xy**) and a nested list of atoms **res_list**.

    :parameter xy: An array with the indices of all *m* core/ligand pairs.
    :type xy: *m*2* |np.ndarray|_ [|np.int64|_]
    :parameter res_list: A list of PLAMS atoms, each nested tuple representing all atoms within
        a given residue.
    :type res_list: |list|_ [|tuple|_ [|plams.Atom|_]]
    :parameter int lig_count: The number of ligand (*n*) in XYn.
    :return:
    """
    ret = {}
    for core, lig in xy.T:
        try:
            ret[res_list[0][core].id].append([at.id for at in res_list[lig]])
        except KeyError:
            ret[res_list[0][core].id] = [[at.id for at in res_list[lig]]]
    for i in ret:
        ret[i] = combinations(ret[i], lig_count)
    return ret
