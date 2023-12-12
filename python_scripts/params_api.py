from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors 
import json
import argparse

app = Flask(__name__)

def get_molecule_params_from_smiles(smiles_string):
   
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        descriptors = Chem.Descriptors.CalcMolDescriptors(mol, missingVal=None, silent=True)
        # descriptors_dict = {desc_name: descriptors[i] for i, desc_name in enumerate(Descriptors._descList)}
       
        return descriptors
    else:
        print("Failed to create a molecule from SMILES.")
        return {}

@app.route('/get_molecule_params', methods=['POST'])


def get_molecule_params():
    try:
        data = request.get_json()
        smiles_string = data.get('smiles')
        params = get_molecule_params_from_smiles(smiles_string)

        if params is not None:
            return jsonify(params), 200
        else:
            return jsonify({'error': 'Failed to create a molecule from SMILES.'}), 400

    except Exception as e:
        print(f"Exception: {e}")
        return jsonify({'error': str(e)}), 400

if __name__ == "__main__":
    app.run(port=5000)

# curl -X POST -H "Content-Type: application/json" -d '{"smiles": "C1=CC(=C(C=C1Cl)C(=O)O)N"}' http://localhost:5000/get_molecule_params
