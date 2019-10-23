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

from collections import abc
from typing import (List, Tuple, Any, Optional)

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import (Molecule, Atom, Settings)

from ..mol_utils import (merge_mol, get_index, round_coords)
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file
from ..workflows.workflow import WorkFlow

__all__ = ['init_qd_construction']

# Aliases for pd.MultiIndex columns
HDF5_INDEX = ('hdf5 index', '')
MOL = ('mol', '')
OPT = ('opt', '')


def init_qd_construction(ligand_df: SettingsDataFrame, core_df: SettingsDataFrame,
                         construct_qd: bool = True) -> SettingsDataFrame:
    """Initialize the quantum dot construction.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of ligands.

    core_df : |CAT.SettingsDataFrame|_
        A dataframe of cores.

    construct_qd : :class:`bool`
        If ``False``, only construct and return the dataframe without filling it with molecules.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    """
    # import pdb; pdb.set_trace()
    qd_df = _get_df(core_df.index, ligand_df.index, ligand_df.settings)
    qd_df[MOL] = None
    qd_df.sort_index(inplace=True)
    if not construct_qd:
        return qd_df

    workflow = WorkFlow.from_template(qd_df, name='qd_attach')
    workflow.keep_files = False

    # Pull from the database; push unoptimized structures
    idx = workflow.from_db(qd_df, inplace=False)

    # Start the ligand optimization
    workflow(construct_mol_series, qd_df, columns=MOL, index=idx,
             core_df=core_df, ligand_df=ligand_df)
    workflow.to_db(qd_df, index=idx)

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    mol_format = qd_df.settings.optional.database.mol_format
    if mol_format and not qd_df.settings.optional.qd.optimize:
        path = workflow.path
        mol_to_file(qd_df.loc[idx, MOL], path, mol_format=mol_format)

    return qd_df


def construct_mol_series(df: SettingsDataFrame, core_df: pd.DataFrame,
                         ligand_df: pd.DataFrame, **kwargs) -> pd.Series:
    """Construct a Series of new quantum dots."""
    def _get_mol(i, j, k, l):
        ij = i, j
        kl = k, l
        return ligand_to_qd(core_df.at[ij, MOL], ligand_df.at[kl, MOL], settings)

    qd_df = df
    settings = qd_df.settings

    mol_list = [_get_mol(i, j, k, l) for i, j, k, l in qd_df.index]
    return pd.Series(mol_list, index=qd_df.index, name=MOL, dtype=object)


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


def ligand_to_qd(core: Molecule, ligand: Molecule, settings: Settings) -> Molecule:
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
    def get_name() -> str:
        core_name = core.properties.name
        anchor = str(qd[-1].properties.pdb_info.ResidueNumber - 1)
        lig_name = ligand.properties.name
        return f'{core_name}__{anchor}_{lig_name}'

    dirname = settings.optional.qd.dirname

    # Define vectors and indices used for rotation and translation the ligands
    vec1 = np.array([-1, 0, 0], dtype=float)  # All ligands are already alligned along the X-axis
    vec2 = np.array(core.get_center_of_mass()) - sanitize_dim_2(core.properties.dummies)
    idx = ligand.get_index(ligand.properties.dummies) - 1
    ligand.properties.dummies.properties.anchor = True

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    lig_array = rot_mol(ligand, vec1, vec2, atoms_other=core.properties.dummies, idx=idx)
    qd = core.copy()
    array_to_qd(ligand, lig_array, mol_out=qd)
    qd.round_coords()

    # Set properties
    qd.properties = Settings({
        'indices': [i for i, at in enumerate(qd, 1) if
                    at.properties.pdb_info.ResidueName == 'COR' or at.properties.anchor],
        'path': dirname,
        'name': get_name(),
        'job_path': [],
        'prm': ligand.properties.prm
    })

    # Print and return
    return qd


def _get_rotmat1(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """Calculate the rotation matrix for rotating **vec1** to **vec2**.

    Returns a unit matrix if the length of one of the vectors is 0.

    Paramaters
    ----------
    vec1 : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.
        Represents the initial orientation.

    vec2 : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.
        Represents the final orientation.

    Returns
    -------
    :math:`n*3*3` |np.ndarray|_
        An array consisting of :math:`n` rotation matrices.

    """
    # Increase the vector dimensionality if required
    _vec1 = np.array(vec1, dtype=float, ndmin=2, copy=False)
    _vec2 = np.array(vec2, dtype=float, ndmin=2, copy=False)

    with np.errstate(divide='raise', invalid='raise'):
        # Return a unit matrix if either vec1 or vec2 is zero
        try:
            u_vec1 = vec1 / np.linalg.norm(_vec1, axis=1)[:, None]
            u_vec2 = vec2 / np.linalg.norm(_vec2, axis=1)[:, None]
        except FloatingPointError:  # i.e. division by zero
            shape = max(len(_vec1), len(_vec2)), 3, 3
            ret = np.zeros(shape, dtype=float)
            ret[:] = np.identity(3)
            return ret

    v1, v2, v3 = np.cross(u_vec1, u_vec2).T
    v0 = np.zeros(max(len(u_vec1), len(u_vec2)))
    M = np.array([[v0, -v3, v2],
                  [v3, v0, -v1],
                  [-v2, v1, v0]]).T

    if len(_vec1) > 1 and len(_vec2) > 1:
        return np.identity(3) + M + ((M@M).T / (1 + u_vec2.T@u_vec1).T[..., None]).T
    elif len(_vec1) == 1:
        return np.identity(3) + M + ((M@M).T / (1 + u_vec2@u_vec1[0].T).T).T
    elif len(_vec2) == 1:
        return np.identity(3) + M + ((M@M).T / (1 + u_vec2[0]@u_vec1.T).T).T
    raise ValueError('vec1 and vec2 expect 1- or 2-dimensional array-like objects; '
                     f'observed shapes: vec1: {np.asarray(vec1).shape} and vec2: '
                     f'{np.asarray(vec2).shape}')


def _get_rotmat2(vec: np.ndarray, step: float = (1/16)) -> np.ndarray:
    r"""Calculate the rotation matrix for rotating m vectors along their axes, each vector
    yielding k = (2 / step) possible rotations.

    Paramaters
    ----------
    vec : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.

    step : float
        The rotation stepsize as fraction of :math:`2*/pi`.

    """
    # Increase the vector dimensionality if required
    _vec = np.array(vec, dtype=float, ndmin=2, copy=False)

    v = _vec / np.linalg.norm(_vec, axis=1)[:, None]
    v1, v2, v3 = v.T
    zero = np.zeros(len(_vec))
    W = np.array([[zero, -v3, v2],
                  [v3, zero, -v1],
                  [-v2, v1, zero]]).T

    step_range = np.pi * np.arange(0.0, 2.0, step)
    a1 = np.sin(step_range)[:, None, None, None]
    a2 = (1 - np.cos(step_range))[:, None, None, None]

    return np.identity(3) + a1 * W + a2 * W@W


def rot_mol(xyz_array: np.ndarray,
            vec1: np.ndarray,
            vec2: np.ndarray,
            idx: int = 0,
            atoms_other: Optional[np.ndarray] = None,
            bond_length: Optional[int] = None,
            step: float = 1/16,
            dist_to_self: bool = True) -> np.ndarray:
    r"""Rotate **xyz_array**.

    Paramaters
    ----------
    xyz_array : :math:`m*n*3` |np.ndarray|_ or :math:`m` |list|_ [|plams.Molecule|_]:
        A 3D array-like object or list of :math:`m` PLAMS molecules consisting of :math:`n` atoms.
        This array represents the to-be rotated molecules.

    vec1 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the initial orientation of **xyz_array**.

    vec2 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the final orientation of **xyz_array**.
        Also used as a rotation axis.

    idx : int:
        An atomic index or sequence of atomic indices.
        Used for translational purposes (see **atoms_other** and **bond_length**).
        Expects 0-based indices.

    atoms_other : :math:`n*3` |np.ndarray|_:
        Optional: A 2D array of Cartesian coordinates.
        If not ``None`` then **xyz_array** will be translated to ensure
        :code:`xyz_array[:, idx] == atoms_other`.

    bond_length : float
        Optional: One or multiple floats representing bond lengths.
        Performs a translation of :code:`xyz_array[:, idx]` along **vec2** by **bond_length**.

    step : float
        The rotation stepsize as fraction of :math:`2*\pi`.

    dist_to_self : bool
        pass

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_:
        The rotated **xyz_array**.

    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz_array = sanitize_dim_3(xyz_array)
    vec1 = sanitize_dim_2(vec1)
    vec2 = sanitize_dim_2(vec2)

    # Define slices
    if xyz_array.ndim == 3 and len(xyz_array) != 1:
        length = None
    else:
        length = max([len(vec1), len(vec2)])
    idx1 = slice(0, None), idx
    idx2 = slice(0, length), slice(int(2 / step)), idx

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
    if dist_to_self is not None or atoms_other is not None:
        a, b, c, d = xyz_array.shape
        ret_array = np.empty((a, c, d), order='F')
        for i, xyz in enumerate(xyz_array):
            dist_array = cdist(xyz.reshape(b*c, d), atoms_other).reshape(b, c, len(atoms_other))
            idx_min = np.nansum(np.exp(-dist_array), axis=(1, 2)).argmin()
            if dist_to_self is not None:
                atoms_other = np.concatenate((atoms_other, xyz[idx_min]))
            ret_array[i] = xyz[idx_min]
        return ret_array
    return xyz_array


def rot_mol_angle(xyz_array: np.ndarray,
                  vec1: np.ndarray,
                  vec2: np.ndarray,
                  idx: int = 0,
                  atoms_other: Optional[np.ndarray] = None,
                  bond_length: Optional[float] = None) -> np.ndarray:
    """Rotate one or more molecules (**xyz_array**) from **vec1** to **vec2**.

    Paramaters
    ----------
    xyz_array : :math:`m*n*3` |np.ndarray|_ or :math:`m` |list|_ [|plams.Molecule|_]:
        A 3D array-like object or list of :math:`m` PLAMS molecules consisting of :math:`n` atoms.
        This array represents the to-be rotated molecules.

    vec1 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the initial orientation of **xyz_array**.

    vec2 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the final orientation of **xyz_array**.

    idx : int:
        An atomic index or sequence of atomic indices.
        Used for translational purposes (see **atoms_other** and **bond_length**).
        Expects 0-based indices.

    atoms_other : :math:`n*3` |np.ndarray|_:
        Optional: A 2D array of Cartesian coordinates.
        If not ``None`` then **xyz_array** will be translated to ensure
        :code:`xyz_array[:, idx] == atoms_other`.

    bond_length : float
        Optional: One or multiple floats representing bond lengths.
        Performs a translation of :code:`xyz_array[:, idx]` along **vec2** by **bond_length**.

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_:
        The rotated **xyz_array**.

    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    _xyz_array = sanitize_dim_3(xyz_array)
    _vec1 = sanitize_dim_2(vec1)
    _vec2 = sanitize_dim_2(vec2)

    # Rotate and translate all n ligands; readjust bond lengths if bond_length is set
    rotmat = _get_rotmat1(_vec1, _vec2)
    ret = _xyz_array@rotmat
    _idx = slice(0, None), idx

    # Translate the the molecules in ret
    if atoms_other is not None:
        atoms_other = sanitize_dim_2(atoms_other)
        ret += (atoms_other - ret[_idx])[:, None, :]
        if bond_length is not None:
            mult = (np.asarray(bond_length) / np.linalg.norm(vec2, axis=1))[:, None]
            ret -= (vec2 * mult)[:, None, :]

    # Return a m*n*3 array
    return ret


def array_to_qd(mol: Molecule, xyz_array: np.ndarray,
                mol_out: Optional[Molecule] = None) -> Optional[List[Molecule]]:
    """Create :math:`n` copies of **mol** and update their Cartesian coordinates with **xyz_array**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A template PLAMS molecule consisting of :math:`n` atoms.

    xyz_array : :math:`m*n*3` |np.ndarray|_:
        A 3D array-like object representing the cartesian coordinates of
        :math:`m` molecules each with :math:`n` atoms.

    mol_out : |plams.Molecule|_
        Optional: If not ``None`` merge the to-be returned molecules with **mol_out**.

    Returns
    -------
    |list|_ [|Molecule|_]
        Optional: if **mol_out** = ``None`` return a list **mol** copies, each with its
        Cartesian coordinates replaced with a set of coordinates from **xyz_array**.

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
    if mol_out is None:
        return mol_list
    else:
        mol_out.merge_mol(mol_list)
        return None


def _is_sequence(item) -> bool:
    return isinstance(item, abc.Sequence)


def _is_atom(item) -> bool:
    return isinstance(item, Atom)


def _is_mol(item) -> bool:
    return isinstance(item, Molecule)


def _is_atom_sequence(item) -> bool:
    return _is_mol(item) or (_is_sequence(item) and _is_atom(item[-1]))


def _is_mol_sequence(item) -> bool:
    return _is_sequence(item) and _is_atom_sequence(item[-1])


def sanitize_dim_2(arg: Any) -> np.ndarray:
    """Convert a PLAMS atom or sequence of :math:`n` PLAMS atoms into a :math:`n*3` array.

    Parameters
    ----------
    arg : object
        The to be parsed object.
        Acceptable types are:

        * A PLAMS :class:`Atom`.
        * A PLAMS :class:`Molecule`.
        * A Sequence consisting of PLAMS :class:`Atom`.
        * An array-like object with a dimensionality smaller than 2 and float-compatible elements.

    Returns
    -------
    :math:`n*3` |np.ndarray|_ [|np.float64|_]
        A 2D array consisting of floats.

    Raises
    ------
    ValueError
        Raised if dimensionality of the to-be returned array is higher than 2
        or the content of **arg** cannot be converted into an array of floats.

    """
    MOL = Molecule()

    if _is_atom(arg):
        return np.array(arg.coords)[None, :]

    elif _is_atom_sequence(arg):
        return MOL.as_array(atom_subset=arg)

    else:
        ret = np.array(arg, ndmin=2, dtype=float, copy=False)
        if ret.ndim > 2:
            raise ValueError(f'Dimensionality of arg ({ret.ndim}) is larger than 2')
        return ret


def sanitize_dim_3(arg: Any,
                   padding: float = np.nan) -> np.ndarray:
    """Convert a Molecule or sequence of :math:`m` molecules into an :math:`m*n*3` array.

    If necessary, the to-be returned array is padded with **padding** .

    Parameters
    ----------
    arg : object
        The to be parsed object.
        Acceptable types are:

        * A PLAMS :class:`Atom`
        * A PLAMS :class:`Molecule`
        * A (nested) Sequences consisting of PLAMS :class:`Atom`
        * A (nested) Sequences consisting of PLAMS :class:`Molecule`
        * An array-like object with a dimensionality smaller than 3 and float-compatible elements.

    padding : float
        A value used for padding the to-be returned array.
        Only relevant if **arg** consists of multiple molecule with different numbers of atoms.

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_ [|np.float64|_]
        A 3D array consisting of floats.

    Raises
    ------
    ValueError
        Raised if dimensionality of the to-be returned array is higher than 3
        or the content of **arg** cannot be converted into an array of floats.

    """
    MOL = Molecule()

    if _is_atom(arg):
        return np.array(arg.coords)[None, None, :]

    elif _is_atom_sequence(arg):
        return MOL.as_array(atom_subset=arg)[None, :]

    elif _is_mol_sequence(arg):
        max_at = max(len(mol) for mol in arg)
        ret = np.full((len(arg), max_at, 3), padding, order='F')
        for i, mol in enumerate(arg):
            ret[i, 0:len(mol)] = MOL.as_array(atom_subset=mol)
        return ret

    else:
        ret = np.array(arg, ndmin=3, dtype=float, copy=False)
        if ret.ndim > 3:
            raise ValueError(f'Dimensionality of arg ({ret.ndim}) is larger than 3')
        return ret
