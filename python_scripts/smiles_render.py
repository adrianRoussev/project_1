from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors 
from rdkit.Chem import Descriptors 

import argparse


def is_reaction_smiles(smiles_string):
    # Check if the input SMILES string represents a reaction
    return ">" in smiles_string

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

def render_reaction_smiles(reaction_smiles, image_file_path):
    reaction = AllChem.ReactionFromSmarts(reaction_smiles)
    if reaction:
        img = Draw.ReactionToImage(reaction)
        try:
            img.save(image_file_path)
            print("Image saved successfully.")
            return True
        except Exception as e:
            print("Error while saving the image:", e)
            return False
    else:
        print("Failed to create a reaction from reaction SMILES.")
        return False

def render_smiles_or_reaction(smiles_string, image_file_path):
    if is_reaction_smiles(smiles_string):
        return render_reaction_smiles(smiles_string, image_file_path)
    else:
        return render_smiles(smiles_string, image_file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Render a molecule or reaction from a SMILES string to an image file.')
    parser.add_argument('smiles', type=str, help='The SMILES string of the molecule or reaction')
    parser.add_argument('image_file', type=str, help='The path to save the image file')
    args = parser.parse_args()

    smiles_string = args.smiles
    image_file_path = args.image_file

    render_smiles_or_reaction(smiles_string, image_file_path)

