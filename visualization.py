# Libraries
import rdkit
from rdkit.Chem import Draw

# Files
import vsepr


def draw_molecule() -> None:
    """Draw 2D representation of molecule using RDKIT and PIL"""
    formula, smiles = vsepr.get_formula_and_smiles()

    molecule = rdkit.Chem.MolFromSmiles(smiles)

    pic = Draw.MolToImage(molecule)

    pic.show()
