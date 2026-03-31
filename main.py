from rdkit import Chem

# ---------------------------------------------------------
# Molecule: Erythromycin
# ---------------------------------------------------------
smiles = "CC[C@H]1[C@@H](C)[C@@H](O)[C@H](OC(=O)[C@H](C)[C@@H](O)[C@H](C)O)[C@@H](C)[C@@H](O)[C@H](C)O[C@H]2[C@H](O)[C@@H](C)O[C@@H](O[C@@H]3[C@H](C)[C@@H](O)[C@H](C)O[C@H]3O)[C@H](O)[C@H]2O[C@@H]1C"
molecule_name = "Erythromycin"

# Create molecule
mol = Chem.MolFromSmiles(smiles)

# Add hydrogens
mol = Chem.AddHs(mol)

# Assign stereochemistry
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

# Find chiral centers
chiral_centers = Chem.FindMolChiralCenters(
    mol,
    includeUnassigned=False,
    useLegacyImplementation=False
)

# Output
print("Molecule:", molecule_name)
print("Total number of chiral centers:", len(chiral_centers))
print()

for atom_index, config in chiral_centers:
    atom = mol.GetAtomWithIdx(atom_index)
    element = atom.GetSymbol()

    print("Atom Index:", atom_index)
    print("Element:", element)
    print("Configuration:", config)
    print()