from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
import json

app = Flask(__name__)

def get_reaction_products(reactants, rxn):
    reaction = AllChem.ReactionFromSmarts(rxn)
    reactant_molecules = [Chem.MolFromSmiles(smiles) for smiles in reactants]
    product_molecules_list = reaction.RunReactants(reactant_molecules)
    
    response = []
    for i, products in enumerate(product_molecules_list):
        reactant_set_result = {"reactant_set": i + 1, "products": []}
        for product in products:
            product_smiles = Chem.MolToSmiles(product)
            reactant_set_result["products"].append({"product_smiles": product_smiles})
        response.append(reactant_set_result)
    
    json_response = json.dumps(response, indent=2)
    return json_response

@app.route('/get_rxn_products', methods=['POST'])
def get_rxn_products():
    try:
        data = request.get_json()
        reactants, rxn = data.get('reactants', []), data.get('rxn', '')
        products = get_reaction_products(reactants, rxn)

        if products is not None:
            return jsonify(products), 200
        else:
            return jsonify({'error': 'Failed to predict products.'}), 400

    except Exception as e:
        print(f"Exception: {e}")
        return jsonify({'error': str(e)}), 400

if __name__ == "__main__":
    app.run(port=6000)


# reagent1 = Chem.MolFromSmiles('NC(=S)N(C)C')
# reagent2 = Chem.MolFromSmiles('[Cl]CC(=O)c1ccccc1')
# rxn_smarts = '[NH2:1][C:2](=[S:3])[NH0:4].[C:5](=[O:6])[C:7][Cl:8]>>[N:4][c:2]1[s:3][c:5][c:7][n:1]1'
# rxn = AllChem.ReactionFromSmarts(rxn_smarts)
# product = rxn.RunReactants((reagent1,reagent2))[0][0]
# Chem.SanitizeMol(product)


# curl -X POST -H "Content-Type: application/json" -d '{
#   "reactants": ["NC(=S)N(C)C", "[Cl]CC(=O)c1ccccc1"],
#   "rxn": "[NH2:1][C:2](=[S:3])[NH0:4].[C:5](=[O:6])[C:7][Cl:8]>>[N:4][c:2]1[s:3][c:5][c:7][n:1]1"
# }' http://localhost:5000/get_rxn_products
