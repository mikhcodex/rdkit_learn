from rdkit import Chem
from icecream import ic
from cheminf_kit import show_molecule

m = Chem.MolFromSmiles('C1=CC2=C(C=C1)C1=CC=CC=C21')
ic(m.GetAtomWithIdx(3).GetIsAromatic())
ic(m.GetAtomWithIdx(6).GetIsAromatic())
ic(m.GetBondBetweenAtoms(3,6).GetIsAromatic())
# show_molecule(m, title="Aromatic example", highlight_atoms=[3,6])
show_molecule(m, title=None, size=(400, 400), kekulize=True, highlight_atoms=None, highlight_bonds=None,
show_atom_indices=False, show_bond_indices=False, block=True)
