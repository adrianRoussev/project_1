# from flask import Flask, request, jsonify
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Draw

# import argparse


# def is_reaction_smiles(smiles_string):
#     # Check if the input SMILES string represents a reaction
#     return ">" in smiles_string

# def render_smiles(smiles_string, image_file_path):
#     mol = Chem.MolFromSmiles(smiles_string)
#     if mol:
#         # Generate 2D coordinates (if not already present)
#         AllChem.Compute2DCoords(mol)
#         img = Draw.MolToImage(mol)
#         try:
#             img.save(image_file_path)
#             print("Image saved successfully.")
#             return True
#         except Exception as e:
#             print("Error while saving the image:", e)
#             return False
#     else:
#         print("Failed to create a molecule from SMILES.")
#         return False

# def render_reaction_smiles(reaction_smiles, image_file_path):
#     reaction = AllChem.ReactionFromSmarts(reaction_smiles)
#     if reaction:
#         img = Draw.ReactionToImage(reaction)
#         try:
#             img.save(image_file_path)
#             print("Image saved successfully.")
#             return True
#         except Exception as e:
#             print("Error while saving the image:", e)
#             return False
#     else:
#         print("Failed to create a reaction from reaction SMILES.")
#         return False

# def render_smiles_or_reaction(smiles_string, image_file_path):
#     if is_reaction_smiles(smiles_string):
#         return render_reaction_smiles(smiles_string, image_file_path)
#     else:
#         return render_smiles(smiles_string, image_file_path)

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='Render a molecule or reaction from a SMILES string to an image file.')
#     parser.add_argument('smiles', type=str, help='The SMILES string of the molecule or reaction')
#     parser.add_argument('image_file', type=str, help='The path to save the image file')
#     args = parser.parse_args()

#     smiles_string = args.smiles
#     image_file_path = args.image_file

#     render_smiles_or_reaction(smiles_string, image_file_path)


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
       
        AllChem.Compute2DCoords(mol)
        try:
            Draw.MolToFile(mol, image_file_path, size=(300, 300), kekulize=True, wedgeBonds=True, imageType='svg', fitImage=False, options=None)
            return True
        except Exception as e:
            print("Error while saving the image:", e)
            return False
    else:
        print("Failed to create a molecule from SMILES.")
        return False
    
def generate_labeled_smarts(reactant_smiles, product_smiles, catalyst_smiles):
    # Convert SMILES to RDKit molecules
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)
    catalyst_mol = Chem.MolFromSmiles(catalyst_smiles)

    # Generate labeled reaction SMARTS
    rxn = AllChem.ReactionFromSmarts(f">> reactant + catalyst >> product")
    rxn.Initialize()
    rxn.AddReactantTemplate(reactant_mol)
    rxn.AddProductTemplate(product_mol)
    rxn.AddAgentTemplate(catalyst_mol)

    # Get the labeled SMARTS
    smarts = AllChem.ReactionToSmarts(rxn)

    return smarts
    
def render_reaction_smiles(reactant_smiles, product_smiles, catalyst_smiles, image_file_path):
    labeled_smarts = generate_labeled_smarts(reactant_smiles, product_smiles, catalyst_smiles)
    reaction = AllChem.ReactionFromSmarts(labeled_smarts, useSmiles = True)
    if reaction:
        d2d = Draw.MolDraw2DCairo(800, 300)
        d2d.DrawReaction(reaction, highlightByReactant=True)
        png = d2d.GetDrawingText()

        try:
            with open(image_file_path, 'wb+').write(png):
                print("Image saved successfully.")
            return True
        except Exception as e:
            print("Error while saving the image:", e)
            return False
    else:
        print("Failed to create a reaction from reaction SMILES.")
        return False


# def render_reaction_smiles(reaction_smiles, image_file_path):
#     reaction = AllChem.ReactionFromSmarts(reaction_smiles)
#     if reaction:
#         try:
#             Draw.ReactionToImage(reaction, image_file_path, subImgSize=(200, 200), useSVG=True, drawOptions=None, returnPNG=False)
#             print("Image saved successfully.")
#             return True
#         except Exception as e:
#             print("Error while saving the image:", e)
#             return False
#     else:
#         print("Failed to create a reaction from reaction SMILES.")
#         return False


def render_smiles_or_reaction(smiles_string, image_file_path):
    if is_reaction_smiles(smiles_string):
        return render_reaction_smiles(reactant_smiles, product_smiles, catalyst_smiles, image_file_path)
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
