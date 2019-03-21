""" A module with misc functions related to manipulating molecules and their geometry. """

__all__ = [
        'merge_mol', 'adf_connectivity', 'fix_h', 'fix_carboxyl',
        'from_mol_other', 'from_rdmol', 'separate_mod'
]

from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.molecule import Molecule
from scm.plams.core.functions import add_to_class
from scm.plams.tools.periodic_table import PeriodicTable
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem
from rdkit.Chem import rdMolTransforms


@add_to_class(Molecule)
def from_mol_other(self, mol, atom_subset=None):
    """ Update the atomic coordinates of *self* with coordinates from another PLAMS molecule.
    Alternatively, update only a subset of atoms.
    Performs an inplace update of **self**.

    :parameter mol: A PLAMS molecule.
    :type mol: |plams.Molecule|_
    :parameter atom_subset: A subset of atoms in **self**.
    :type atom_subset: |None|_ or |list|_ [|plams.Atom|_]
    """
    atom_subset = atom_subset or self.atoms
    for at1, at2 in zip(atom_subset, mol):
        at1.coords = at2.coords


@add_to_class(Molecule)
def from_rdmol(self, rdmol, atom_subset=None):
    """ Update the atomic coordinates of *self* with coordinates from an RDKit molecule.
    Alternatively, update only a subset of atoms.
    Performs an inplace update of **self**.

    :parameter rdmol: An RDKit molecule.
    :type rdmol: |rdkit.Chem.Mol|_
    :parameter atom_subset: A subset of atoms in **self**.
    :type atom_subset: |None|_ or |list|_ [|plams.Atom|_]
    """
    atom_subset = atom_subset or self.atoms
    conf = rdmol.GetConformer()
    for at1, at2 in zip(atom_subset, rdmol.GetAtoms()):
        pos = conf.GetAtomPosition(at2.GetIdx())
        at1.coords = (pos.x, pos.y, pos.z)


def to_atnum(item):
    """ Turn an atomic symbol into an atomic number.

    :parameter item: An atomic symbol or number.
    :type item: |int|_ or |str|_
    :return: An atomic number.
    :rtype: |int|_
    """
    if isinstance(item, str):
        return PeriodicTable.get_atomic_number(item)
    return item


def to_symbol(item):
    """ Turn an atomic number into an atomic symbol.

    :parameter item: An atomic symbol or number.
    :type item: |int|_ or |str|_
    :return: An atomic symbol.
    :rtype: |str|_
    """
    if isinstance(item, int):
        return PeriodicTable.get_symbol(item)
    return item


@add_to_class(Atom)
def get_atom_index(self):
    """ Return the index of an atom (numbering starts with 1).

    :return: An atomic index.
    :rtype: |int|_.
    """
    return self.mol.atoms.index(self) + 1


@add_to_class(Bond)
def get_bond_index(self):
    """ Return a tuple of two atomic indices defining a bond (numbering starts with 1).

    :return: A tuple of 2 atomic indices defining a bond.
    :rtype: 2 |tuple|_ [|int|_].
    """
    return self.atom1.get_atom_index(), self.atom2.get_atom_index()


@add_to_class(Molecule)
def merge_mol(self, mol_list):
    """ Merge two or more molecules into a single molecule.
    No new copies of atoms/bonds are created, all atoms/bonds are moved from mol_list to plams_mol.
    Performs an inplace update of **self**.

    :parameter mol_list: A molecule or list of molecules.
    :type mol_list: |plams.Molecule|_ or |list|_ [|plams.Molecule|_].
    """
    if isinstance(mol_list, Molecule):
        mol_list = [mol_list]

    for mol in mol_list:
        for atom in mol.atoms:
            atom.mol = self
        for bond in mol.bonds:
            bond.mol = self
        self.properties.soft_update(mol.properties)
        self.atoms += mol.atoms
        self.bonds += mol.bonds


@add_to_class(Molecule)
def separate_mod(self):
    """ Modified PLAMS function: seperates a molecule instead of a copy of a molecule.
    Separate the molecule into connected components.
    returns is a list of new |Molecule| objects (all atoms and bonds are disjoint with
        the original molecule).
    Each element of this list is identical to one connected component of the base molecule.
    A connected component is a subset of atoms such that there exists a path
        (along one or more bonds) between any two atoms.

    :return: A list of molecules with atoms and bonds from **self**.
    :rtype: |list|_ [|plams.Molecule|_]
    """
    frags = []
    for at in self:
        at._visited = False

    def dfs(v, mol):
        v._visited = True
        v.mol = mol
        for e in v.bonds:
            e.mol = mol
            u = e.other_end(v)
            if not u._visited:
                dfs(u, mol)

    for src in self.atoms:
        if not src._visited:
            m = Molecule()
            dfs(src, m)
            frags.append(m)

    for at in self.atoms:
        del at._visited
        at.mol.atoms.append(at)
    for b in self.bonds:
        b.mol.bonds.append(b)

    return frags


def adf_connectivity(plams_mol):
    """ Create an AMS-compatible connectivity list.

    :parameter plams_mol: A PLAMS molecule.
    :type plams_mol: |plams.Molecule|_
    :return: An ADF-compatible connectivity list of *n* bonds.
    :rtype: *n* |list|_ [|str|_].
    """
    # Create list of indices of all aromatic bonds
    rdmol = molkit.to_rdmol(plams_mol)
    aromatic = [bond.GetIsAromatic() for bond in rdmol.GetBonds()]

    # Create a list of bond orders; aromatic bonds get a bond order of 1.5
    plams_mol.set_atoms_id()
    bond_orders = [bond.order for bond in plams_mol.bonds]
    for i, ar in enumerate(aromatic):
        if ar:
            bond_orders[i] = 1.5
    bonds = [str(bond.atom1.id) + ' ' + str(bond.atom2.id) + ' ' + str(order) for
             bond, order in zip(plams_mol.bonds, bond_orders)]
    plams_mol.unset_atoms_id()

    return bonds


def fix_carboxyl(plams_mol):
    """ Resets carboxylate OCO angles if it is smaller than 60 degrees.
    Performs an inplace update of **plams_mol**.

    :parameter plams_mol: A PLAMS molecule.
    :type plams_mol: |plams.Molecule|_
    """
    rdmol = molkit.to_rdmol(plams_mol)
    carboxylate = Chem.MolFromSmarts('[O-]C(C)=O')
    matches = rdmol.GetSubstructMatches(carboxylate)

    if matches:
        get_angle = rdMolTransforms.GetAngleDeg
        set_angle = rdMolTransforms.SetAngleDeg
        for idx in matches:
            if get_angle(rdmol.GetConformer(), idx[3], idx[1], idx[0]) < 60:
                set_angle(rdmol.GetConformer(), idx[2], idx[1], idx[3], 180.0)
                set_angle(rdmol.GetConformer(), idx[0], idx[1], idx[3], 120.0)
        plams_mol.from_rdmol(rdmol)


def fix_h(plams_mol):
    """ If a C=C-H angle is smaller than 20.0 degrees, set it back to 120.0 degrees.
    Performs an inplace update of **plams_mol**.

    :parameter plams_mol: A PLAMS molecule.
    :type plams_mol: |plams.Molecule|_
    """
    H_list = [atom for atom in plams_mol if atom.atnum == 1 and 2.0 in
              [bond.order for bond in plams_mol.neighbors(atom)[0].bonds]]

    rdmol = molkit.to_rdmol(plams_mol)
    idx = plams_mol.atoms.index
    set_angle = rdMolTransforms.SetAngleDeg
    get_angle = rdMolTransforms.GetAngleDeg

    update = []
    for atom in H_list:
        at1 = atom
        at2 = plams_mol.neighbors(at1)[0]
        at3 = [atom for atom in plams_mol.neighbors(at2) if atom != at1]
        if get_angle(rdmol.GetConformer(), idx(at3[0]), idx(at2), idx(at1)) <= 20.0:
            set_angle(rdmol.GetConformer(), idx(at3[0]), idx(at2), idx(at1), 120.0)
            update.append(True)
        elif get_angle(rdmol.GetConformer(), idx(at3[1]), idx(at2), idx(at1)) <= 20.0:
            set_angle(rdmol.GetConformer(), idx(at3[1]), idx(at2), idx(at1), 120.0)
            update.append(True)

    if update:
        plams_mol.from_rdmol(rdmol)
