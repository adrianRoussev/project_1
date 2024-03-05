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

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=1500)
