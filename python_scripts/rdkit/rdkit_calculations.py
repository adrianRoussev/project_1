
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
from rdkit.Chem import FunctionalGroups
import rdkit.Chem.Fragments as Fragments
from rdkit.Chem import rdFMCS
import pandas as pd



class ReactionHandler:

    def __init__(self):
        csv_file_path = "fragment_data_sheet2.csv"
        self.df = pd.read_csv(csv_file_path)
        

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
    
    def get_compounds_from_smarts(self, rxn_smarts, compound_type):
        try:
            reaction = AllChem.ReactionFromSmarts(rxn_smarts)
            if reaction is not None:
                if compound_type == "reactants":
                    compounds = reaction.GetReactants()
                elif compound_type == "catalysts":
                    compounds = reaction.GetAgents()
                elif compound_type == "products":
                    compounds = reaction.GetProducts()
                else:
                    return None
                
                return [Chem.MolToSmiles(compound) for compound in compounds]
                
                # if mol_type == "smiles":
                #     
                # if mol_type == "smarts":
                #     return [Chem.MolToSmarts(compound) for compound in compounds]
                # if mol_type == "mol":
                #     return compounds
    
            else:
                return None
        except:
            return None
    
    @staticmethod
    def is_smarts_or_smiles(string):
        try:
            if Chem.MolFromSmarts(string) is not None:
                return "SMARTS"
            elif Chem.MolFromSmiles(string) is not None:
                return "SMILES"
            else:
                return None
        except:
            return None
        
    def calculate_molecule_similarity(self, mol_1, mol_2):
        # Compute 2D coordinates and SSSR for both molecules
        AllChem.Compute2DCoords(mol_1)
        Chem.GetSSSR(mol_1)

        AllChem.Compute2DCoords(mol_2)
        Chem.GetSSSR(mol_2)

        # Check for substructure match
        substructure_match = mol_1.HasSubstructMatch(mol_2).float

        return float(substructure_match)
    
    def calculate_score(self, reactant, product, target):
        # params = rdFMCS.MCSParameters()
        # params.AtomCompare = rdFMCS.AtomCompare.CompareElements
        # params.BondCompare = rdFMCS.BondCompare.CompareAny
        # params.RingMatchesRingOnly = False

        # Find MCS for the first check
        mcs = rdFMCS.FindMCS([reactant, product])
        matched_atoms = mcs.smartsString
        matched_atoms_count = mcs.numAtoms
        if matched_atoms_count > 1:
            score = (matched_atoms_count/target.GetNumAtoms() )
        else:
            score = 0
        
        return score

    
    def find_reacting_atoms(self, rxn_smarts):
        
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        reactants = [reactant for reactant in reaction.GetReactants()]
        products = [product for product in reaction.GetProducts()]

        reactant_hash = {}

        for reactant in reactants: 
            reactant_smiles = Chem.MolToSmarts(reactant)
            list = []
            for product in products:
                product_smiles = Chem.MolToSmarts(product)
                reactant_similarity = self.calculate_score(reactant, product, reactant)
                product_similarity = self.calculate_score(reactant, product, product)
                # product_hash = {}
                # product_hash[product_smiles] = similarity

                list.append({"smiles": product_smiles, "similarity": reactant_similarity, "product_similarity": product_similarity})
                

            list = sorted(list, key=lambda x: x.get('similarity', 0), reverse=True)
            reactant_hash[reactant_smiles] = list
        
        print(reactant_hash)
        return reactant_hash
        
    




            
        #         reactant_hash[(reactant, product)] = similarity

        # sorted_hash = sorted(reactant_hash.items(), key=lambda x: x[1], reverse=True)
        
        # print(sorted_hash)
        # reacting_atoms_list = []
        
        # for reactant_product, score in sorted_hash:
                
        #         reactant = reactant_product[0]
        #         product = reactant_product[1]

                  
        #         reacting_atoms = set()
    
        #         reacting_atoms.update(
        #         set(range(reactant.GetNumAtoms())).symmetric_difference(
        #         set(range(product.GetNumAtoms()))))
       
        #         reacting_atoms_list.append(reacting_atoms)

        #         print(reacting_atoms, Chem.MolToSmiles(product), Chem.MolToSmiles(reactant),score)
      
    
    def get_mol_functional_groups(self, mol):
        all_functions = [func for func in dir(Fragments) if callable(getattr(Fragments, func))]

        responses = {}
        value = {}
        for func_name in all_functions:
            smarts = self.get_functional_group_smarts(func_name)
            func = getattr(Fragments, func_name)
            try:
                response = func(mol)
                responses[func_name] = response
                if response != 0:
                    value[smarts] = response 
            except Exception as e:
                responses[func_name] = str(e)

        smarts_list = []

        for smarts in value.keys():
            if smarts is not None: 
                smarts_list.append(smarts)
        
        return smarts_list



    def get_functional_group_smarts(self, function_name):
        filtered_rows = self.df[self.df['Code'] == function_name]
    
        if not filtered_rows.empty:
            smarts_value = filtered_rows.iloc[0]['SMARTS']
            description = filtered_rows.iloc[0]['Description']

            smarts_and_description = (smarts_value, description)
            return smarts_and_description
        


    def find_reaction_functional_groups(self, smarts):
        
        rxn_smarts = self.get_standardized_smarts(smarts)

        reactants = []
        for each in self.get_compounds_from_smarts(rxn_smarts, "reactants"):
            reactant = Chem.MolFromSmiles(each)
            reactants.append(reactant)

        products = []
        for each in self.get_compounds_from_smarts(rxn_smarts, "products"):
            product = Chem.MolFromSmiles(each)
            products.append(product)

        reactant_fgs = []
        product_fgs = []
        fg_dict = {}
        for reactant in reactants:
            fgs = self.get_mol_functional_groups(reactant)
            fgs.append(Chem.MolToSmarts(reactant))
            reactant_fgs.append(fgs)
           
        for product in products:
            fgs = self.get_mol_functional_groups(product)
            fgs.append(Chem.MolToSmarts(product))
            product_fgs.append(fgs)
    
        fg_dict['reactants'] = reactant_fgs
        fg_dict['products'] = product_fgs

        # print(fg_dict)
        return fg_dict

        
        
    def find_changing_functional_groups(self, smarts):   

        fg_dict = self.find_reaction_functional_groups(smarts)
        reactant_fgs = fg_dict['reactants'] 
        product_fgs = fg_dict['products']
        
        changing_groups = []
        new_groups = []
        reaction_fgs = {}

        for reactant_fg in reactant_fgs:
            if reactant_fg[0] not in [each[0] for each in product_fgs]:
                changing_groups.append(reactant_fg)

        for product_fg in product_fgs:
            if product_fg[0] not in [each[0] for each in reactant_fgs]:
                new_groups.append(product_fg)
                
        reaction_fgs['changing reactant groups'] = changing_groups
        reaction_fgs['new product groups'] = new_groups

        return reaction_fgs
    
    def get_base_rxn(self, smarts, reagents=False ):
        reaction_fgs = self.find_changing_functional_groups(smarts)

        base_reactant_patterns = []

        for list in reaction_fgs['changing reactant groups']:
            base_reactant = []
            for item in list:
                if isinstance(item, tuple):
                    base_reactant.append(item[0])
            base_reactant_patterns.append(base_reactant)
                
        base_product_patterns = []

        for list in reaction_fgs['new product groups']:
            base_product = []
            for item in list:
                if isinstance(item, tuple):
                    base_product.append(item[0])
            base_product_patterns.append(base_product)

        reactants = [each[-1] for each in reaction_fgs['changing reactant groups']]
        products = [each[-1] for each in reaction_fgs['new product groups']]

        print({'reactant_patterns': base_reactant_patterns, 'reactant': reactants })
        print( {'product_patterns': base_product_patterns, 'product': products})

        reactant_count = range(len(reactants))
        reactant_mcs = []

        for each in reactant_count:
            for pattern in base_reactant_patterns[each]:
                
                mcs = rdFMCS.FindMCS([Chem.MolFromSmarts(reactants[each]), Chem.MolFromSmarts(pattern)])
                if mcs.smartsString != '':
                    reactant_mcs.append(Chem.MolFromSmarts(mcs.smartsString))
        
        product_count = range(len(products))
        product_mcs = []

        for each in product_count:
            for pattern in base_product_patterns[each]:
                
                mcs = rdFMCS.FindMCS([Chem.MolFromSmarts(products[each]), Chem.MolFromSmarts(pattern)])
                if mcs.smartsString != '':
                    product_mcs.append(Chem.MolFromSmarts(mcs.smartsString))
        

        reaction = Chem.rdChemReactions.ChemicalReaction()

        for reactant in reactant_mcs:
            if reactant is not None:
                reaction.AddReactantTemplate(reactant)

        for product in product_mcs:
            if product is not None:
                reaction.AddProductTemplate(product)
        
        if reagents:
            original_reaction = AllChem.ReactionFromSmarts(smarts) 
            reagent_list = [Chem.MolToSmiles(agent) for agent in original_reaction.GetAgents()]

            print({"REAGENTS": reagent_list})
            for reagent in reagent_list:
                reagent_mol = Chem.MolFromSmiles(reagent)
                reaction.AddAgentTemplate(reagent_mol)
            
        reaction_smarts = self.get_standardized_smarts(AllChem.ReactionToSmarts(reaction))
        print("Reaction SMARTS:", reaction_smarts)
        return reaction_smarts

    def run_base_rxn(self, smarts, reagents=False):

        if reagents:
            base_rxn_smarts = self.get_base_rxn(smarts, True)

        else:
            base_rxn_smarts = self.get_base_rxn(smarts)

        print({'base_rxn': base_rxn_smarts,
                'rxn': smarts})
        
        base_rxn = AllChem.ReactionFromSmarts(base_rxn_smarts)
        reaction = AllChem.ReactionFromSmarts(smarts)

        reactants = [Chem.MolToSmiles(reactant) for reactant in reaction.GetReactants()]
        reactant_smiles  = [Chem.MolFromSmiles(smiles) for smiles in reactants]
        reactant_tuple = tuple(reactant_smiles)
        
        print(reactant_tuple)
        result_products = base_rxn.RunReactants(reactant_tuple)

        if result_products:
            print({'results': result_products})
            results = []
            pair = []
            for each in result_products:
                for product in each:
                    pair.append(Chem.MolToSmarts(product))
            results.append(pair)

            print({"results": results})
            return results
        else:
            print({"results": "reaction did not run"})
            return None
        

        



    




        

        






        



        
        # mol = Chem.MolFromSmiles(mol_string)

        # hierarchy = FunctionalGroups.BuildFuncGroupHierarchy()

        # # all_nodes = []

        # # for level in hierarchy:
        # #     for node in level:
        # #         all_nodes.append(node)

        # # for node in all_nodes:
        # #     print(node)

        # functional_groups = []
        # for level, group in enumerate(hierarchy):
        #     pattern = group.rxnSmarts 
            
            # for each in pattern:
            #     print(each.smarts)
            # match = mol.GetSubstructMatch(Chem.MolFromSmarts(pattern))
            # print(match)
        #     functional_groups.append({'Name': group['Name'], 'Atoms': match})

        # print("Functional Groups:")
        # for fg in functional_groups:
        #     print(f"{fg['Name']}: {fg['Atoms']}")

# Example usage:
# reaction_handler = ReactionHandler()
rxn_smarts ='[#6](-[#6])-[#8].[#6]-[#6](=[#8])-[#1]>>[#6]-[#6](=[#8])-[#8]-[#6]-[#6]'


# rxn_dict = reaction_handler.get_mols_from_rxn(rxn_smarts)
# reaction = reaction_handler.get_rxn_from_mols(rxn_dict)

# reactant1 = 'CC=O'
# reactant2 = 'CC(=O)CC'
# product1 = 'CCOC(=O)C'

mol_string = 'C1=CC(=CC=C1C(=O)O)Cl'

reaction_handler = ReactionHandler()
# reaction_handler.get_mol_functional_groups(mol_string)
# reaction_handler.find_changing_functional_groups(rxn_smarts)
# reaction_handler.get_base_rxn(rxn_smarts, True)
reaction_handler.run_base_rxn(rxn_smarts)

# reaction_handler.find_reacting_atoms(rxn_smarts)


# reactants = [reactant1, reactant2]
# reactant = reactant1
# products = [product1, product2]
# catalysts = []
# rxn_dict = {
#         'reactants': reactants,
#         'products': products,
#         'catalysts': catalysts
#     }
# 
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
    
    