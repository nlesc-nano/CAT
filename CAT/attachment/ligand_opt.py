""" A module designed for optimizing the geometry of ligands. """

__all__ = ['init_ligand_opt']

import itertools
import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.core.errors import MoleculeError
from scm.plams.core.settings import Settings
from scm.plams.core.functions import add_to_class
from scm.plams.tools.units import Units
from scm.plams.recipes.global_minimum import global_minimum_scan_rdkit
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit.Chem import AllChem

from .ligand_attach import (rot_mol_angle, sanitize_dim_2)
from ..data_handling.CAT_database import Database
from ..utils import get_time
from ..mol_utils import (to_symbol, fix_carboxyl, get_bond_index,
                         from_mol_other, from_rdmol, separate_mod)
from ..data_handling.database import (mol_from_database, mol_to_database)


def init_ligand_opt(ligand_df, arg):
    """ Initialize the ligand optimization procedure.
    Performs an inplace update of **ligand_df**.

    :parameter ligand_df: A dataframe of valid ligands.
    :type ligand_df: |pd.DataFrame|_ (columns: |str|_, index=|int|_, values=|plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    """
    database = Database(path=arg.optional.database.dirname)
    overwrite = 'ligand' in arg.optional.database.overwrite:

    # Searches for matches between the input ligand and the database; imports the structure
    if 'ligand' in arg.optional.database.read:
        database.from_csv(ligand_df, database='ligand')
        for i, mol in zip(ligand_df['hdf5 index'], ligand_df['mol']):
            if i >= 0:
                print(get_time() + mol.properties.name + '\t has been pulled from the database')
        print('')

    # Optimize all new ligands
    if arg.optional.ligand.optimize:
        # Identify the to be optimized ligands
        if overwrite:
            idx = ligand_df.index
            message = '\t has been (re-)optimized'
        else:
            idx = ligand_df['hdf5 index'] < 0
            message = '\t has been optimized'

        # Optimize the ligands
        for i, ligand in ligand_df['mol'][idx].iteritems():
            mol_list = split_mol(ligand)
            for mol in mol_list:
                mol.set_dihed(180.0)
            ligand_tmp = recombine_mol(mol_list)
            fix_carboxyl(ligand_tmp)

            # Print messages
            print(get_time() + ligand.properties.name + message)
        print('')

    # Write newly optimized structures to the database
    if 'ligand' in arg.optional.database.write and arg.optional.ligand.optimize:
        recipe = Settings()
        recipe.settings = {'name': '1', 'key': 'RDKit', 'value': 'UFF'}
        overwrite = 'ligand' in arg.optional
        database.update_csv(ligand_df, columns=['formula', 'hdf5 index'],
                            job_recipe=recipe, database='ligand', overwrite=overwrite)


@add_to_class(Molecule)
def split_bond(self, bond, atom_type='H', bond_length=1.1):
    """ Delete a bond and cap the resulting fragments.
    A link to the two atoms previously defining the bond & the two capping atoms is stored under
        self.properties.mark in a list of 4-tuples.
    Performs in inplace update of **self**.

    :parameter bond: A PLAMS bond.
    :type:  |plams.Bond|_
    :parameter atom_type: The atomic symbol or number of the two to be created capping atoms.
    :type atom_type: |str|_ or |int|_
    :parameter float bond_length: The length of the two new bonds in angstrom.
    """
    atom_type = to_symbol(atom_type)
    at1, at2 = bond.atom1, bond.atom2
    at3, at4 = Atom(symbol=atom_type, coords=at1.coords), Atom(symbol=atom_type, coords=at2.coords)
    self.add_atom(at3, adjacent=[at2])
    self.add_atom(at4, adjacent=[at1])
    self.bonds[-1].resize(at1, bond_length)
    self.bonds[-2].resize(at2, bond_length)
    if self.properties.mark:
        self.properties.mark.append((at1, at4, at2, at3))
    else:
        self.properties.mark = [(at1, at4, at2, at3)]
    self.delete_bond(bond)


@add_to_class(Molecule)
def neighbors_mod(self, atom, exclude=1):
    """ A modified PLAMS function: Allows the exlucison of specific elements from the return list.
    Return a list of neighbors of **atom** within the molecule. Atoms with
    **atom** has to belong to the molecule. Returned list follows the same order as the
    **atom.bond** attribute.

    :parameter atom: The plams atom whose nieghbours will be returned.
    :type atom: |plams.Atom|_
    :parameter int exclude: Exclude all neighbours with a specific atomic number.
    :return: A list of all neighbours of **atom**.
    :rtype: |list|_ [|plams.Atom|_].
    """
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum != exclude]


def split_mol(plams_mol):
    """ Split a molecule into multiple smaller fragments,
    one fragment for every branch within **plams_mol**.

    :parameter plams_mol: The input molecule with the properties.dummies attribute.
    :type plams_mol: |plams.Molecule|_
    :return: A list of one or more plams molecules.
    :rtype: |list|_ [|plams.Molecule|_]
    """
    # Temporary remove hydrogen atoms
    h_atoms = []
    h_bonds = []
    for atom in reversed(plams_mol.atoms):
        if atom.atnum == 1:
            h_atoms.append(atom)
            h_bonds.append(atom.bonds[0])
            plams_mol.delete_atom(atom)

    # Remove undesired bonds
    bond_list = [bond for bond in plams_mol.bonds if not plams_mol.in_ring(bond.atom1) and not
                 plams_mol.in_ring(bond.atom2)]

    # Remove even more undesired bonds
    for bond in reversed(bond_list):
        n1, n2 = plams_mol.neighbors_mod(bond.atom1), plams_mol.neighbors_mod(bond.atom2)
        if not (len(n1) >= 3 and len(n2) >= 2) and not (len(n1) >= 2 and len(n2) >= 3):
            bond_list.remove(bond)

    # Add the hydrogen atoms and bonds back to the molecule
    for atom, bond in zip(h_atoms, h_bonds):
        plams_mol.add_atom(atom)
        plams_mol.add_bond(bond)

    atom_list = list(itertools.chain.from_iterable((bond.atom1, bond.atom2) for bond in bond_list))
    atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
    atom_dict = {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}

    # Fragment the molecule such that the functional group is on the largest fragment
    for at in atom_dict:
        for i in atom_dict[at][2:]:
            len_atom = [plams_mol.get_frag_size(bond, plams_mol.properties.dummies) for
                        bond in atom_dict[at]]
            idx = len_atom.index(max(len_atom))
            bond = atom_dict[at][idx]
            plams_mol.split_bond(bond)
            atom_dict[at].remove(bond)

    # Copy the properties attribute to all fragment molecules
    properties = plams_mol.properties
    mol_list = plams_mol.separate_mod()
    for mol in mol_list:
        mol.properties = properties

    return mol_list


@add_to_class(Molecule)
def get_frag_size(self, bond, atom):
    """ Return the size of a moleculair fragment containing **atom** if **self** was split into two
    molecules by the breaking of **bond**.

    :parameter bond: A PLAMS bond.
    :type bond: |plams.Bond|_
    :parameter atom: A PLAMS atom. The size of the fragment containg this atom will be returned.
    :type atom: |plams.Atom|_
    :return: The number of atoms in the fragment containing **atom**.
    :rtype: |int|_.
    """
    if bond not in self.bonds:
        error = 'get_frag_size: The argument bond should be of type plams.Bond and be part'
        error += ' of the Molecule'
        raise MoleculeError(error)
    elif atom not in self.atoms:
        error = 'get_frag_size: The argument atom should be of type plams.Atom and be part'
        error += ' of the Molecule'
        raise MoleculeError(error)

    for at in self:
        at._visited = False

    def dfs(at1, ar=np.zeros(2), atom=atom):
        at1._visited = True
        ar[0] += 1
        if at1 is atom:
            ar[1] = 1
        for bond in at1.bonds:
            at2 = bond.other_end(at1)
            if not at2._visited:
                ar += dfs(at2)
        return ar

    bond.atom1._visited = bond.atom2._visited = True
    size1, has_atom1 = dfs(bond.atom1)
    size2, has_atom2 = dfs(bond.atom2)
    for at in self.atoms:
        del at._visited

    if has_atom1:
        return size1
    return size2


def recombine_mol(mol_list):
    """ Recombine a list of molecules into a single molecule.
    A list of 4-tuples of plams.Atoms will be read from mol_list[0].properties.mark.
    A bond will be created between tuple[0] & tuple[2]; tuple[1] and tuple[3] will be deleted.

    :parameter mol_list: A list of on or more plams molecules with the properties.mark atribute.
    :type: |list|_ [|plams.Molecule|_]
    :return: The (re-)merged PLAMS molecule.
    :rtype: |plams.Molecule|_.
    """
    if len(mol_list) == 1:
        return mol_list[0]
    tup_list = mol_list[0].properties.mark
    if not tup_list:
        error = 'No PLAMS atoms specified in mol_list[0].properties.mark, aborting recombine_mol()'
        raise IndexError(error)

    for tup in tup_list:
        # Allign mol1 & mol2
        mol1, mol2 = tup[0].mol, tup[2].mol
        vec1 = sanitize_dim_2(tup[3]) - sanitize_dim_2(tup[2])
        vec2 = sanitize_dim_2(tup[0]) - sanitize_dim_2(tup[1])
        idx = tup[2].get_atom_index() - 1
        mol_array = rot_mol_angle(mol2, vec1, vec2, atoms_other=tup[0], idx=idx, bond_length=1.5)
        mol2.from_array(mol_array)

        # Merge mol1 & mol2
        mol1.merge_mol(mol2)
        mol1.delete_atom(tup[1])
        mol1.delete_atom(tup[3])
        mol1.add_bond(tup[0], tup[2])
        bond_tup = mol1.bonds[-1].get_bond_index()
        mol1.from_mol_other(global_minimum_scan_rdkit(mol1, bond_tup))

    del mol1.properties.mark
    return mol1


def get_dihed(atoms, unit='degree'):
    """ Returns the dihedral angle defined by four atoms.

    :parameter atoms: An iterable consisting of 4 PLAMS atoms
    :type atoms: 4 |tuple|_ [|plams.atoms|_]
    :parameter str unit: The output unit.
    :return: A dihedral angle in **unit**.
    :rtype: |float|_.
    """
    vec1 = -np.array(atoms[0].vector_to(atoms[1]))
    vec2 = np.array(atoms[1].vector_to(atoms[2]))
    vec3 = np.array(atoms[2].vector_to(atoms[3]))

    v1v2, v2v3 = np.cross(vec1, vec2), np.cross(vec3, vec2)
    v1v2_v2v3 = np.cross(v1v2, v2v3)
    v2_norm_v2 = vec2 / np.linalg.norm(vec2)
    epsilon = np.arctan2(v1v2_v2v3@v2_norm_v2, v1v2@v2v3)

    return Units.convert(epsilon, 'radian', unit)


@add_to_class(Molecule)
def set_dihed(self, angle, opt=True, unit='degree'):
    """ Change a dihedral angle into a specific value.
    Performs an inplace update of **self**.

    :parameter float angle: The desired dihedral angle.
    :parameter bool opt: Whether or not the dihedral adjustment should be followed up by an
        RDKit UFF optimization.
    :parameter str unit: The input unit.
    """
    angle = Units.convert(angle, unit, 'degree')
    bond_list = [bond for bond in self.bonds if bond.atom1.atnum != 1 and bond.atom2.atnum != 1
                 and bond.order == 1 and not self.in_ring(bond)]

    for bond in bond_list:
        n1, n2 = self.neighbors_mod(bond.atom1), self.neighbors_mod(bond.atom2)
        n1 = [atom for atom in n1 if atom != bond.atom2]
        n2 = [atom for atom in n2 if atom != bond.atom1]
        if len(n1) > 1:
            n1 = [atom for atom in n1 if len(self.neighbors_mod(atom)) > 1]
        if len(n2) > 1:
            n2 = [atom for atom in n2 if len(self.neighbors_mod(atom)) > 1]
        if n1 and n2:
            dihed = get_dihed((n1[0], bond.atom1, bond.atom2, n2[0]))
            if self.properties.dummies not in bond:
                self.rotate_bond(bond, bond.atom1, angle - dihed, unit='degree')
            else:
                self.rotate_bond(bond, bond.atom1, -dihed, unit='degree')

    if opt:
        rdmol = molkit.to_rdmol(self)
        AllChem.UFFGetMoleculeForceField(rdmol).Minimize()
        self.from_rdmol(rdmol)
