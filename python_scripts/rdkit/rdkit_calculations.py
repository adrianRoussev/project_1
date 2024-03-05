
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

class ReactionHandler:
    def get_mols_from_rxn(self, rxn_smarts):
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        Chem.rdChemReactions.SanitizeRxn(reaction)
        
        reactants = [Chem.MolToSmiles(reactant) for reactant in reaction.GetReactants()]
        products = [Chem.MolToSmiles(product) for product in reaction.GetProducts()]
        catalysts = [Chem.MolToSmiles(catalyst) for catalyst in reaction.GetAgents()]

        rxn_dict = {
            'reactants': reactants,
            'products': products,
            'catalysts': catalysts
        }

        return rxn_dict

    def get_indexed_smarts(self, smiles_string):
        mol = Chem.MolFromSmiles(smiles_string)
       
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        
        smarts_with_indices = Chem.MolToSmarts(mol)
        return smarts_with_indices

    def get_standardized_smarts(self, rxn_smarts):
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        Chem.rdChemReactions.SanitizeRxn(reaction)
        
        reactants = [reactant for reactant in reaction.GetReactants()]
        products = [product for product in reaction.GetProducts()]
        catalysts = [catalyst for catalyst in reaction.GetAgents()]

        reaction = rdChemReactions.ChemicalReaction()    
        total_atom_count = 0

        for mol in reactants:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1 + total_atom_count) 
            total_atom_count += mol.GetNumAtoms() 
            reaction.AddReactantTemplate(mol)

        for mol in catalysts:
            reaction.AddAgentTemplate(mol)

        total_atom_count = 0  

        for mol in products:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1 + total_atom_count)  
            total_atom_count += mol.GetNumAtoms() 
            reaction.AddProductTemplate(mol)

        return rdChemReactions.ReactionToSmarts(reaction)

# Example usage:
# reaction_handler = ReactionHandler()
# rxn_smarts = "[OH:1][CH3:2].[CH3:3][C:4](=[O:5])[O:6][H:7]>>[O:1]=[C:2]([OH:3])[CH3:4].[O:5]=[C:6]=[O:7]"
# rxn_dict = reaction_handler.get_mols_from_rxn(rxn_smarts)
# reaction = reaction_handler.get_rxn_from_mols(rxn_dict)

# reactant1 = 'CC=O'
# reactant2 = 'CC(=O)CC'
# product1 = 'CCOC(=O)C'
# product2 = 'CCOCC(=O)CC'

# reactants = [reactant1, reactant2]
# reactant = reactant1
# products = [product1, product2]
# catalysts = []
# rxn_dict = {
#         'reactants': reactants,
#         'products': products,
#         'catalysts': catalysts
#     }
# reaction_handler = ReactionHandler()
# rxn = reaction_handler.get_rxn_from_mols(rxn_dict)
# print(rdChemReactions.ReactionToSmarts(rxn))


# curl -X POST \
#   http://localhost:1500/get_rxn_from_mols \
#   -H 'Content-Type: application/json' \
#   -d '{
#         "rxn_dict": {
#             "reactants": ["CC=O", "CC(=O)CC"],
#             "products": ["CCOC(=O)C", "CCOCC(=O)CC"],
#             "catalysts": []
#         }
#     }'
    
    