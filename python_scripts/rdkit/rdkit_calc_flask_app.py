from flask import Flask, request, jsonify
from rdkit_calculations import ReactionHandler
import json

app = Flask(__name__)
reaction_handler = ReactionHandler()

@app.route('/get_mols_from_rxn', methods=['POST'])
def get_mols_from_rxn():
    data = request.get_json()
    rxn_smarts = data['rxn_smarts']
    rxn_dict = reaction_handler.get_mols_from_rxn(rxn_smarts)
    rxn_dict_json = json.dumps(rxn_dict)  # Serialize the dictionary to JSON
    return jsonify(rxn_dict_json)

@app.route('/get_indexed_smarts', methods=['POST'])
def get_indexed_smarts():
    data = request.get_json()
    smiles_string = data['smiles_string']
    smarts_with_indices = reaction_handler.get_indexed_smarts(smiles_string)
    return jsonify({'smarts_with_indices': smarts_with_indices})

@app.route('/get_standardized_smarts', methods=['POST'])
def get_standardized_smarts():
    data = request.get_json()
    rxn_smarts = data['rxn_smarts']
    standardized_smarts = reaction_handler.get_standardized_smarts(rxn_smarts)
    return jsonify({'standardized_smarts': standardized_smarts})

@app.route('/get_compounds_from_smarts', methods=['POST'])
def get_compounds_from_smarts():
    data = request.get_json()
    rxn_smarts = data['rxn_smarts']
    compound_type = data['compound_type']
    compounds = reaction_handler.get_compounds_from_smarts(rxn_smarts, compound_type)
    return jsonify({'compounds': compounds})

@app.route('/get_mol_functional_groups', methods=['POST'])
def get_mol_functional_groups():
    data = request.get_json()
    smiles_string = data['smiles_string']
    mol = Chem.MolFromSmiles(smiles_string)
    smarts_list = reaction_handler.get_mol_functional_groups(mol)
    return jsonify({'functional_groups': smarts_list})

@app.route('/find_reaction_functional_groups', methods=['POST'])
def find_reaction_functional_groups():
    data = request.get_json()
    smarts = data['smarts']
    functional_groups = reaction_handler.find_reaction_functional_groups(smarts)
    return jsonify({'functional_groups': functional_groups})

@app.route('/find_changing_functional_groups', methods=['POST'])
def find_changing_functional_groups():
    data = request.get_json()
    smarts = data['smarts']
    changing_groups = reaction_handler.find_changing_functional_groups(smarts)
    return jsonify({'changing_groups': changing_groups})

@app.route('/get_base_rxn', methods=['POST'])
def get_base_rxn():
    data = request.get_json()
    smarts = data['smarts']
    reagents = data.get('reagents', False)
    base_rxn_smarts = reaction_handler.get_base_rxn(smarts, reagents)
    return jsonify({'base_rxn_smarts': base_rxn_smarts})

@app.route('/run_base_rxn', methods=['POST'])
def run_base_rxn():
    data = request.get_json()
    smarts = data['smarts']
    reagents = data.get('reagents', False)
    results = reaction_handler.run_base_rxn(smarts, reagents)
    return jsonify({'results': results})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=1500)
