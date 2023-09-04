from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import argparse

def render_smiles(smiles_string, image_file_path):
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        # Generate 2D coordinates (if not already present)
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol)
        try:
            img.save(image_file_path)
            print("Image saved successfully.")
            return True
        except Exception as e:
            print("Error while saving the image:", e)
            return False
    else:
        print("Failed to create a molecule from SMILES.")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Render a molecule from a SMILES string to an image file.')
    parser.add_argument('smiles', type=str, help='The SMILES string of the molecule')
    parser.add_argument('image_file', type=str, help='The path to save the image file')
    args = parser.parse_args()

    smiles_string = args.smiles
    image_file_path = args.image_file

    render_smiles(smiles_string, image_file_path)
