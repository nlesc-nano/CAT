"""
CAT.attachment.ligand_attach
============================

A module designed for attaching ligands to cores.

Index
-----
.. currentmodule:: CAT.attachment.ligand_attach
.. autosummary::
    init_qd_construction
    construct_mol_series
    _read_database
    _get_indices
    _get_df
    ligand_to_qd
    _get_rotmat1
    _get_rotmat2
    rot_mol
    rot_mol_angle
    array_to_qd
    sanitize_dim_2
    sanitize_dim_3

API
---
.. autofunction:: init_qd_construction
.. autofunction:: construct_mol_series
.. autofunction:: _read_database
.. autofunction:: _get_indices
.. autofunction:: _get_df
.. autofunction:: ligand_to_qd
.. autofunction:: _get_rotmat1
.. autofunction:: _get_rotmat2
.. autofunction:: rot_mol
.. autofunction:: rot_mol_angle
.. autofunction:: array_to_qd
.. autofunction:: sanitize_dim_2
.. autofunction:: sanitize_dim_3

"""

from typing import (List, Tuple)

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams.mol.molecule import Molecule
from scm.plams.core.settings import Settings

from ..settings_dataframe import SettingsDataFrame
from ..logger import logger
from ..mol_utils import (merge_mol, get_index)

try:
    from dataCAT import (Database, mol_to_file)
    DATA_CAT = True
except ImportError:
    DATA_CAT = False

__all__ = ['init_qd_construction']

# Aliases for pd.MultiIndex columns
HDF5_INDEX = ('hdf5 index', '')
MOL = ('mol', '')
OPT = ('opt', '')


def init_qd_construction(ligand_df: SettingsDataFrame,
                         core_df: SettingsDataFrame) -> SettingsDataFrame:
    """Initialize the quantum dot construction.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of ligands.

    core_df : |CAT.SettingsDataFrame|_
        A dataframe of cores.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    """
    # Extract arguments
    settings = ligand_df.settings.optional
    overwrite = DATA_CAT and 'qd' in settings.database.overwrite
    write = DATA_CAT and 'qd' in settings.database.write
    read = DATA_CAT and 'qd' in settings.database.read
    qd_path = settings.qd.dirname
    db_path = settings.database.dirname
    mol_format = settings.database.mol_format

    # Attempt to pull structures from the database
    qd_df = _get_df(core_df.index, ligand_df.index, ligand_df.settings)
    qd_df.sort_index(inplace=True)
    if read:
        mol_series1 = _read_database(qd_df, ligand_df, core_df)

    # Identify and create the to be constructed quantum dots
    mol_series2 = construct_mol_series(qd_df, core_df, ligand_df)

    # Update the *mol* column in qd_df with 1 or 2 series of quantum dots
    try:
        qd_df[MOL] = mol_series1.append(mol_series2)
    except NameError:
        qd_df[MOL] = mol_series2

    # Export the resulting geometries back to the database
    if write:
        data = Database(db_path, **settings.database.mongodb)
        data.update_csv(qd_df, columns=[HDF5_INDEX], database='QD_no_opt')
        mol_to_file(qd_df[MOL], qd_path, overwrite, mol_format)
    return qd_df


def construct_mol_series(qd_df: SettingsDataFrame,
                         core_df: pd.DataFrame,
                         ligand_df: pd.DataFrame) -> pd.Series:
    """Construct a Series of new quantum dots"""
    def _get_mol(i, j, k, l):
        ij = i, j
        kl = k, l
        return ligand_to_qd(core_df.at[ij, MOL], ligand_df.at[kl, MOL], settings)

    settings = qd_df.settings
    idx = qd_df[HDF5_INDEX] < 0

    mol_list = [_get_mol(i, j, k, l) for i, j, k, l in qd_df.index[idx]]
    return pd.Series(mol_list, index=qd_df.index[idx], name=MOL, dtype=object)


def _read_database(qd_df: SettingsDataFrame,
                   ligand_df: SettingsDataFrame,
                   core_df: SettingsDataFrame) -> pd.Series:
    """Read quantum dots from the database and set their properties.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of ligands.

    core_df : |CAT.SettingsDataFrame|_
        A dataframe of cores.

    Returns
    -------
    |pd.Series|_ [|plams.Molecule|_]
        A Series of quantum dots pulled from the database.

    """
    def get_name():
        """Construct the name of a quantum dot."""
        core = core_df.at[(i[0:2]), MOL].properties.name
        res = mol[-1].properties.pdb_info.ResidueNumber - 1
        lig = ligand_df.at[(i[2:4]), MOL].properties.name
        return '{}__{:d}_{}'.format(core, res, lig)

    # Extract arguments
    settings = qd_df.settings.optional
    path = settings.database.dirname
    data = Database(path, **settings.database.mongodb)

    # Extract molecules from the database and set their properties
    # If possible extract optimized structures; supplement with unoptimized structures if required
    mol_series_opt = data.from_csv(qd_df, database='QD', inplace=False)
    mol_series_no_opt = data.from_csv(qd_df, database='QD_no_opt', inplace=False)
    slice_ = mol_series_no_opt.index.isin(mol_series_opt.index)
    mol_series = mol_series_opt.append(mol_series_no_opt[~slice_])

    # Update Molecule.properties
    logger.info('Pulling quantum dots from database')
    for i, mol in mol_series.iteritems():
        mol.properties = Settings({
            'indices': _get_indices(mol, i),
            'path': path,
            'job_path': [],
            'name': get_name()
        })
        logger.info(f'{mol.properties.name} has been pulled from the database')
    return mol_series


def _get_indices(mol: Molecule,
                 index: Tuple[str, str, str, str]) -> List[int]:
    """Return a list with the indices of all atoms in the core plus ligand anchor atoms.

    Ligand anchor atoms are furthermore marked with the properties.anchor attribute.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    index : |tuple|_ [|str|_]
        A tuple of 4 strings.

    Returns
    -------
    |list|_ [|int|_]
        A list of atomic indices.

    """
    # Collect the indices of the atoms in the core
    ret = []
    for i, at in enumerate(mol, 1):
        if at.properties.pdb_info.ResidueName == 'COR':
            ret.append(i)
        else:
            break

    # Extract the index (within the ligand) of the ligand anchor atom
    index = index[3]
    for j, _ in enumerate(index):
        try:
            k = int(index[j:]) - 1
            break
        except ValueError:
            pass
    k += i - 1

    # Append and return
    ref_name = mol[k+1].properties.pdb_info.Name
    for i, at in enumerate(mol.atoms[k:], k+1):
        if at.properties.pdb_info.Name == ref_name:
            at.properties.anchor = True
            ret.append(i)
    return ret


def _get_df(core_index: pd.MultiIndex,
            ligand_index: pd.MultiIndex,
            settings: Settings) -> SettingsDataFrame:
    """Create and return a new quantum dot dataframe.

    Parameters
    ----------
    core_index : |pd.MultiIndex|_
        A multiindex of the cores.

    ligand_index : |pd.MultiIndex|_
        A multiindex of the ligands.

    settings : |plams.Settings|_
        A Settings intance extracted from the ligand or core dataframe.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        An empty dataframe (*i.e.* filled with ``-1`` and ``False``) of quantum dots.

    """
    # Create the index
    idx_tups = [(i, j, k, l) for i, j in core_index for k, l in ligand_index]
    index = pd.MultiIndex.from_tuples(
        idx_tups, names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )

    # Create the collumns
    column_tups = [HDF5_INDEX, OPT]
    columns = pd.MultiIndex.from_tuples(column_tups, names=['index', 'sub index'])

    # Create and return the quantum dot dataframe
    data = {HDF5_INDEX: -1, OPT: False}
    return SettingsDataFrame(data, index=index, columns=columns, settings=settings)


def ligand_to_qd(core: Molecule,
                 ligand: Molecule,
                 settings: Settings) -> Molecule:
    """Function that handles quantum dot (qd, *i.e.* core + all ligands) operations.

    Combine the core and ligands and assign properties to the quantom dot.

    Parameters
    ----------
    core : |plams.Molecule|_
        A core molecule.

    ligand : |plams.Molecule|_
        A ligand molecule.

    settings : |plams.Settings|_
        A settings object containing all (optional) arguments.

    Returns
    -------
    |plams.Molecule|_
        A quantum dot consisting of a core molecule and *n* ligands

    """
    def get_name():
        ret = core.properties.name + '__'
        ret += str(qd[-1].properties.pdb_info.ResidueNumber - 1) + '_' + ligand.properties.name
        return ret

    dirname = settings.optional.qd.dirname

    # Define vectors and indices used for rotation and translation the ligands
    vec1 = sanitize_dim_2(ligand.properties.dummies) - np.array(ligand.get_center_of_mass())
    vec2 = np.array(core.get_center_of_mass()) - sanitize_dim_2(core.properties.dummies)
    idx = ligand.get_index(ligand.properties.dummies) - 1
    ligand.properties.dummies.properties.anchor = True

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    lig_array = rot_mol(ligand, vec1, vec2, atoms_other=core.properties.dummies, idx=idx)
    qd = core.copy()
    array_to_qd(ligand, lig_array, mol_other=qd)

    # Set properties
    qd.properties = Settings({
        'indices': [i for i, at in enumerate(qd, 1) if
                    at.properties.pdb_info.ResidueName == 'COR' or at.properties.anchor],
        'path': dirname,
        'name': get_name(),
        'job_path': []
    })

    # Print and return
    return qd


def _get_rotmat1(vec1, vec2):
    """ Calculate the rotation matrix for rotating m1 vectors in vec1 to m2 vectors in vec2.
    The following values of m1 & m2 are acceptable:
        m1 = m2
        m1 = 1, m2 > 1
        m1 > 1, m2 = 1 """
    with np.errstate(divide='raise', invalid='raise'):
        # Return a unit matrix if either vec1 or vec2 is zero
        try:
            a = vec1 / np.linalg.norm(vec1)
            b = vec2 / np.linalg.norm(vec2, axis=1)[:, None]
        except FloatingPointError:
            shape = max(len(vec1), len(vec2)), 3, 3
            ret = np.zeros(shape, dtype=float)
            ret[:] = np.identity(3)
            return ret

    v1, v2, v3 = np.cross(a, b).T
    zero = np.zeros(len(b))
    M = np.array([[zero, -v3, v2],
                  [v3, zero, -v1],
                  [-v2, v1, zero]]).T
    return np.identity(3) + M + ((M@M).T / (1 + b@a.T).T).T


def _get_rotmat2(vec, step=(1/16)):
    """ Calculate the rotation matrix for rotating m vectors along their axes, each vector
    yielding k = (2 / step) possible rotations. """
    v = vec / np.linalg.norm(vec, axis=1)[:, None]
    v1, v2, v3 = v.T
    zero = np.zeros(len(vec))
    W = np.array([[zero, -v3, v2],
                  [v3, zero, -v1],
                  [-v2, v1, zero]]).T

    step_range = np.pi * np.arange(0.0, 2.0, step)
    a1 = np.sin(step_range)[:, None, None, None]
    a2 = (1 - np.cos(step_range))[:, None, None, None]
    return np.identity(3) + a1 * W + a2 * W@W


def rot_mol(xyz_array, vec1, vec2, idx=0, atoms_other=None, bond_length=False, step=1/16,
            dist_to_self=True):
    """  Define a m*3*3 rotation matrix using vec1 and vec2.
    Depending on the dimension of xyz_array, vec1 & vec2 the following operations can be conducted:
        Rotate 1 molecule by 1 rotation matrix    > 1 rotated molecule
        Rotate 1 molecule by m rotation matrices  > m copies of the mol at different rotations
        Rotate m molecules by 1 rotation matrix   > m molecules rotated in an identical fashion
        Rotate m molecules by m rotation matrices > m molecules, each rotated differently

    Numpy arrays and (nested) iterable consisting of PLAMS atoms can be used interchangeably for
    xyz_array and mol_other; Atomic coordinates will be extracted if necessary and cast into an
    appropriately shaped array.

    :parameter xyz_array: An 3d array or list of PLAMS to be rotated PLAMS molecules consisting of
        *m* molecules with up to *n* atoms.
    :type xyz_array: *m*n*3* |np.ndarray|_ [|np.float64|_] or *m* |list|_ [|plams.Molecule|_]
    :parameter vec1: One or multiple vectors representing the initial orientation of **xyz_array**.
    :type vec1: *m*3* or *3* |np.ndarray|_ [|np.float64|_]
    :parameter vec1: One or multiple vectors representing the final orientation of **xyz_array**.
    :type vec2: *m*3* or *3* |np.ndarray|_ [|np.float64|_]
    :parameter idx: An atomic index or iterable consisting of multiple atomic indices. Translate
        **xyz_array**[:, **idx**] to **atoms_other**.
    :type idx: |int|_ or *m* |np.ndarray|_ [|np.int64|_]
    :parameter atoms_other: A list of PLAMS atoms or a 2d array. Translate
        **xyz_array**[:, **idx**] to **atoms_other**.
    :type atoms_other: |None|_, *m* |list|_ [|plams.Atom|_] or *m*3* |np.ndarray|_ [|np.float64|_]
    :parameter bond_length: The distance(s) between **idx** and **atoms_other**.
    :type bond_length: |float|_ or *m* |np.ndarray|_ [|np.float64|_]
    :return: An array with *m* rotated molecules with up to *n* atoms.
    :rtype: *m*n*3* or *n*3* |np.ndarray|_ [|np.float64|_].
    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz_array = sanitize_dim_3(xyz_array)
    vec1 = sanitize_dim_2(vec1)
    vec2 = sanitize_dim_2(vec2)

    # Define slices
    if xyz_array.ndim == 3 and len(xyz_array) != 1:
        length = len(xyz_array)
    else:
        length = max([len(vec1), len(vec2)])
    idx1 = np.arange(len(xyz_array)), idx
    idx2 = np.arange(length), slice(int(2 / step)), idx

    # Translate xyz_array[idx] to the origin and rotate
    xyz_array -= xyz_array[idx1]
    rotmat = _get_rotmat1(vec1, vec2)
    xyz_array = xyz_array@rotmat

    # Create all k possible rotations of all m ligands
    rotmat = _get_rotmat2(vec2, step=step)
    xyz_array = np.swapaxes(xyz_array@rotmat, 0, 1)

    # Translate the the molecules in xyz_array
    if atoms_other is not None:
        atoms_other = sanitize_dim_2(atoms_other)
        xyz_array += (atoms_other[:, None, :] - xyz_array[idx2])[..., None, :]
        if bond_length:
            mult = (np.asarray(bond_length) / np.linalg.norm(vec2, axis=1))[:, None]
            xyz_array -= (vec2 * mult)[:, None, :]

    # Returns the conformation of each molecule that maximizes the inter-moleculair distance
    # Or return all conformations if dist_to_self = False and atoms_other = None
    if dist_to_self or atoms_other is not None:
        a, b, c, d = xyz_array.shape
        ret_array = np.empty((a, c, d), order='F')
        for i, xyz in enumerate(xyz_array):
            dist_array = cdist(xyz.reshape(b*c, d), atoms_other).reshape(b, c, len(atoms_other))
            idx_min = np.nansum(np.exp(-dist_array), axis=(1, 2)).argmin()
            if dist_to_self:
                atoms_other = np.concatenate((atoms_other, xyz[idx_min]))
            ret_array[i] = xyz[idx_min]
        return ret_array
    return xyz_array


def rot_mol_angle(xyz_array, vec1, vec2, idx=0, atoms_other=None, bond_length=False):
    """
    Define a m*3*3 rotation matrix using vec1 and vec2.
    Depending on the dimension of xyz_array, vec1 & vec2 the following operations can be conducted:
        Rotate 1 molecule by 1 rotation matrix    > 1 rotated molecule
        Rotate 1 molecule by m rotation matrices  > m copies of the mol at different rotations
        Rotate m molecules by 1 rotation matrix   > m molecules rotated in an identical fashion
        Rotate m molecules by m rotation matrices > m molecules, each rotated differently
    Numpy arrays and (nested) iterable consisting of PLAMS atoms can be used interchangeably for
    xyz_array and mol_other; Atomic coordinates will be extracted if necessary and cast into an
    appropriately shaped array.
    xyz_array <np.ndarray>: A m*n*3 array, a n*3 numpy array or a (nested) iterable consisting of
        the to be rotated PLAMS atoms.
    vec1 & vec2 <np.ndarray>: Two n*3 and/or 3 arrays representing one or more vectors.
        vec1 defines the initial vector(s) and vec2 defines the final vector(s) adopted by the to
        be rotated molecule(s).
    atoms_other <None> or <plams.Molecule>: None or a n*3 array, a 3 numpy array, a PLAMS atom or
        iterable consisting of PLAMS atoms. All molecules will be translated to ensure that the
        atom with the index idx adopts the same position as mol_other.
    idx <int>: An atomic index or iterable consisting of multiple atomic indices.
    bond_length <float>: A float or iterable consisting of floats. After translation from
        xyz_array[:, idx] to mol_other, xyz_array will be further translated along vec2 by
        a user-specified bond length.
    return <np.ndarray>: A m*n*3 array or n*3 array representing xyz coordinates of m rotated
        molecules, each consisting of n atoms.
    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz_array = sanitize_dim_3(xyz_array)
    vec1 = sanitize_dim_2(vec1)
    vec2 = sanitize_dim_2(vec2)

    # Rotate and translate all n ligands; readjust bond lengths if bond_length is set
    rotmat = _get_rotmat1(vec1, vec2)
    xyz_array = xyz_array@rotmat
    idx = np.arange(xyz_array.shape[0]), idx

    # Translate the the molecules in xyz_array
    if atoms_other is not None:
        atoms_other = sanitize_dim_2(atoms_other)
        xyz_array += (atoms_other - xyz_array[idx])[:, None, :]
        if bond_length:
            mult = (np.asarray(bond_length) / np.linalg.norm(vec2, axis=1))[:, None]
            xyz_array -= (vec2 * mult)[:, None, :]

    # Return a n*3 or m*n*3 array
    if xyz_array.shape[0] == 1:
        return xyz_array[0]
    return xyz_array


def array_to_qd(mol, xyz_array, mol_other=False):
    """
    Takes a template molecule with n atoms and create m copies of this molecule, updating its
        cartesian coordinates based on a m*n*3 array.
    Has the option to directly combine all copied molecules with a second molecule (mol_other)
        instead of returning the rotated molecules.

    mol <plams.Molecule>: A template PLAMS molecule consisting of n atoms.
    xyz_array <np.ndarray>: A m*n*3 or n*2 numpy array representing the cartesian coordinates of
        m molecules each with n atoms.
    mol_other <plams.Molecule> or False: Add all atoms and bonds from the to be returned molecule
        list to this molecule, instead of returning the molecule list.
    return <list> or None: Returns a list of rotated molecules or concatenates aforementioned
        list and add its atoms and bonds to mol_other, returning nothing.
    """
    mol_list = []
    xyz_array = sanitize_dim_3(xyz_array)
    for i, xyz in enumerate(xyz_array, 2):
        mol_cp = mol.copy()
        mol_cp.from_array(xyz)
        for at in mol_cp:
            at.properties.pdb_info.ResidueNumber = i
        mol_list.append(mol_cp)

    # return a list of molecules or merge the list with an existing molecule
    if not mol_other:
        return mol_list
    mol_other.merge_mol(mol_list)


def sanitize_dim_2(arg):
    """
    Convert a PLAMS atom or iterable consisting of n PLAMS atoms into a n*3 array.
    In addition, 3 arrays are converted into n*3 arrays.
    """
    if not isinstance(arg, np.ndarray):
        try:
            return np.array(arg.coords)[None, :]
        except AttributeError:
            dummy = Molecule()
            return dummy.as_array(atom_subset=arg)
    else:
        if len(arg.shape) == 1:
            return arg[None, :]
        return arg


def sanitize_dim_3(arg, padding=np.nan):
    """
    Convert an iterable consisting of n PLAMS atoms or a nested iterable consisting of m*(≤n)
    PLAMS atoms into a m*n*3 array, padding the array n is not constant. In addition, n*3 arrays
    are converted into m*n*3 arrays.
    """
    if not isinstance(arg, np.ndarray):
        dummy = Molecule()
        try:
            return dummy.as_array(atom_subset=arg)[None, :, :]
        except AttributeError:
            max_at = max(len(mol) for mol in arg)
            ret = np.empty((len(arg), max_at, 3), order='F')
            ret[:] = padding
            for i, mol in enumerate(arg):
                ret[i, 0:len(mol)] = dummy.as_array(atom_subset=mol)
            return ret
    else:
        if len(arg.shape) == 2:
            return arg[None, :, :]
        return arg
