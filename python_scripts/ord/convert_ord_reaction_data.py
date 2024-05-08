from ord_schema.proto import dataset_pb2
from rdkit import Chem
from rdkit.Chem import AllChem


class ReactionSimilarity:
    def __init__(self, file, input_smiles, cutoff_score, similarity_type, params):
        self.file = file
        self.input_smiles = input_smiles
        self.cutoff_score = cutoff_score
        self.similarity_type = similarity_type
        self.params = params

    def parse_data(self):
        with open(self.file, 'rb') as f:
            data = f.read()

        dataset = dataset_pb2.Dataset()
        dataset.ParseFromString(data)

        values_list = []

        for reaction in dataset.reactions:
            reaction_identifiers = reaction.identifiers
            for identifier in reaction_identifiers:
                value = identifier.value
                values_list.append(value)
        return values_list

    def get_reaction_ids(self):
        with open(self.file, 'rb') as f:
            data = f.read()

        dataset = dataset_pb2.Dataset()
        dataset.ParseFromString(data)

        reaction_id_list = []

        for reaction in dataset.reactions:
            reaction_ids = reaction.reaction_id
            reaction_id_list.append(reaction_ids)

        return reaction_id_list

    def calculate_molecule_similarity(self, mol_1, mol_2):
        # Compute 2D coordinates and SSSR for both molecules
        AllChem.Compute2DCoords(mol_1)
        Chem.GetSSSR(mol_1)

        AllChem.Compute2DCoords(mol_2)
        Chem.GetSSSR(mol_2)

        # Check for substructure match
        substructure_match = mol_1.HasSubstructMatch(mol_2)

        return float(substructure_match)

    def calculate_set_similarity(self, set_1, set_2):
        similarities = []
        for mol_1 in set_1:
            for mol_2 in set_2:
                similarity = self.calculate_molecule_similarity(mol_1, mol_2)
                similarities.append(similarity)

        total_similarity = sum(similarities)

        if len(set_1) > len(set_2):
            count = len(set_1)
        elif len(set_2) > len(set_1):
            count = len(set_2)
        else:
            count = len(set_1)

        score = total_similarity / count

        return score

    def calculate_reaction_similarity(self, smiles1, smiles2):
        reaction_1 = AllChem.ReactionFromSmarts(smiles1)
        reaction_2 = AllChem.ReactionFromSmarts(smiles2)

        reactants_1 = list(reaction_1.GetReactants())
        products_1 = list(reaction_1.GetProducts())

        reactants_2 = list(reaction_2.GetReactants())
        products_2 = list(reaction_2.GetProducts())

        for mol in reactants_1 + products_1 + reactants_2 + products_2:
            AllChem.Compute2DCoords(mol)
            Chem.GetSSSR(mol)

        total_similarity = []
        reactants_similarity = self.calculate_set_similarity(reactants_1, reactants_2)

        products_similarity = self.calculate_set_similarity(products_1, products_2)
        product_count = len(products_1) + len(products_2) / 2
        reactant_count = len(reactants_1) + len(reactants_2) / 2
        product_ratio = product_count / (product_count + reactant_count)
        reactant_ratio = reactant_count / (product_count + reactant_count)
        total_similarity.append(reactants_similarity * reactant_ratio)
        total_similarity.append(products_similarity * product_ratio)
        similarity_score = sum(total_similarity)

        similarity_dict = {}
        similarity_dict['reactant similarity'] = reactants_similarity
        similarity_dict['product similarity'] = products_similarity
        similarity_dict['total similarity'] = similarity_score

        return similarity_dict

    def calculate_similarity_from_file(self):
        values_list = self.parse_data()
        results_list = []

        # Perform similarity search for each value in values_list
        for value in values_list:
            result = self.calculate_reaction_similarity(self.input_smiles, value)
            results_list.append(result)
        return results_list

    def get_top_similar_reactions_with_ids(self):
        values_list = self.parse_data()
        results_list = []

        # Perform similarity search for each value in values_list
        for reaction_id, value in zip(self.get_reaction_ids(), values_list):
            result = self.calculate_reaction_similarity(self.input_smiles, value)

            if result[self.similarity_type] >= self.cutoff_score:
                results_list.append((reaction_id, result))

        return results_list

    def get_top_similar_reactions_with_param(self):
        values_list = self.parse_data()
        params_result_list = []

        # Perform similarity search for each value in values_list
        for parameter, value in zip(self.get_reaction_param(), values_list):
            result = self.calculate_reaction_similarity(self.input_smiles, value)

            if result[self.similarity_type] >= self.cutoff_score:
                params_result_list.append((parameter, result))

        params_result_list = sorted(params_result_list, key=lambda x: x[1][self.similarity_type], reverse=True)

        return params_result_list

    def get_reaction_param(self):
        with open(self.file, 'rb') as f:
            data = f.read()

        dataset = dataset_pb2.Dataset()
        dataset.ParseFromString(data)

        reaction_param_list = []
        reaction_param_dict = {}
        for reaction in dataset.reactions:
            for param in self.params:
                reaction_params = getattr(reaction, param)
                reaction_param_dict[param] = reaction_params
                reaction_param_list.append(reaction_param_dict)

        return reaction_param_list

# Example usage
# input_smiles = '[CH3:1][NH:2][CH3:3].Cl[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1>>[CH3:1][N:2]([CH3:3])[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1'

# file = '/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb'
# similarity_type = 'total similarity'
# cutoff_score = 0.5
# params = ['identifiers', 'reaction_id']

# reaction_data = ReactionSimilarity(file, input_smiles, cutoff_score, similarity_type, params)



# result_top_similar_reactions_with_param = reaction_data.get_top_similar_reactions_with_param()
# print(result_top_similar_reactions_with_param)

# curl -X POST -H "Content-Type: application/json" -d '{"file": "/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb","input_smiles": "[CH3:1][NH:2][CH3:3].Cl[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1>>[CH3:1][N:2]([CH3:3])[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1","cutoff_score": 0.5,"similarity_type": "total similarity","params": ["identifiers", "reaction_id"],"function_name": "calculate_similarity_from_file"}' http://localhost:9000/reaction-similarity


# from ord_schema.proto import dataset_pb2
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import DataStructs
# from rdkit.DataStructs import FingerprintSimilarity

# def parse_data(file_path):
#     with open(file_path, 'rb') as f:
#         data = f.read()

#     dataset = dataset_pb2.Dataset()
#     dataset.ParseFromString(data)

#     values_list = []
    
#     for reaction in dataset.reactions:
#         reaction_identifiers = reaction.identifiers
#         for identifier in reaction_identifiers:
#             value = identifier.value
#             values_list.append(value)
#     return values_list

# def get_reaction_ids(file_path):
#     with open(file_path, 'rb') as f:
#         data = f.read()

#     dataset = dataset_pb2.Dataset()
#     dataset.ParseFromString(data)

#     reaction_id_list = []
    
#     for reaction in dataset.reactions:
#         reaction_ids = reaction.reaction_id
            
#         reaction_id_list.append(reaction_ids)
    
#     return reaction_id_list

# def get_reaction_param(file_path, params):
#     with open(file_path, 'rb') as f:
#         data = f.read()

#     dataset = dataset_pb2.Dataset()
#     dataset.ParseFromString(data)
    
#     reaction_param_list = []
#     reaction_param_dict= {}
#     for reaction in dataset.reactions:
#         for param in params:
#             reaction_params = getattr(reaction, param)
#             reaction_param_dict[param] = reaction_params
#             reaction_param_list.append(reaction_param_dict)
            
    
#     return reaction_param_list

# def calculate_molecule_similarity(mol_1, mol_2):
#     # Compute 2D coordinates and SSSR for both molecules
#     AllChem.Compute2DCoords(mol_1)
#     Chem.GetSSSR(mol_1)
    
#     AllChem.Compute2DCoords(mol_2)
#     Chem.GetSSSR(mol_2)

#     # Check for substructure match
#     substructure_match = mol_1.HasSubstructMatch(mol_2)

#     return float(substructure_match)

# def calculate_set_similarity(set_1, set_2):
#     similarities = []
#     for mol_1 in set_1:
#         for mol_2 in set_2:
#             similarity = calculate_molecule_similarity(mol_1, mol_2)
#             similarities.append(similarity)
            
#             total_similarity = sum(similarities)

#     if len(set_1) > len(set_2):
#         count = len(set_1)
#     elif len(set_2)>len(set_1):
#         count= len(set_2)
#     else: 
#         count = len(set_1)
  
#     score = total_similarity / count
   
#     return score
            

# def calculate_reaction_similarity(smiles1, smiles2):
#     reaction_1 = AllChem.ReactionFromSmarts(smiles1)
#     reaction_2 = AllChem.ReactionFromSmarts(smiles2)

#     reactants_1 = list(reaction_1.GetReactants())
#     products_1 = list(reaction_1.GetProducts())

#     reactants_2 = list(reaction_2.GetReactants())
#     products_2 = list(reaction_2.GetProducts())
    
#     for mol in reactants_1 + products_1 + reactants_2 + products_2:
#         AllChem.Compute2DCoords(mol)
#         Chem.GetSSSR(mol)

#     total_similarity = []
#     reactants_similarity =  calculate_set_similarity(reactants_1, reactants_2)
   
#     products_similarity=  calculate_set_similarity(products_1, products_2)
#     product_count = len(products_1) +len(products_2)/2
#     reactant_count = len(reactants_1) +len(reactants_2)/2
#     product_ratio = product_count/(product_count + reactant_count)
#     reactant_ratio = reactant_count/(product_count + reactant_count)
#     total_similarity.append(reactants_similarity*reactant_ratio)
#     total_similarity.append(products_similarity*product_ratio)
#     similarity_score = sum(total_similarity)
    
#     similarity_dict = {}
#     similarity_dict['reactant similarity'] = reactants_similarity
#     similarity_dict['product similarity'] = products_similarity
#     similarity_dict['total similarity'] = similarity_score
    
#     return similarity_dict
#     # reactant_similarity = calculate_similarity(smarts_1, smarts_2)

    
#     # reactant_similarity_list = []
#     # reactant_similarity_list.append(reactant_similarity)

# def calculate_similarity_from_file(file, input_smiles):
#     values_list = parse_data(file)
#     results_list = []
#     # Perform similarity search for each value in values_list
#     for value in values_list:
#         result=  calculate_reaction_similarity(input_smiles, value)
#         results_list.append(result)
#     return results_list
   

# # def get_top_similar_reactions(file, input_smiles, cutoff_score, similarity_type):
# #     values_list = parse_data(file)
# #     results_list = []
# #     # Perform similarity search for each value in values_list
# #     for value in values_list:
# #         result=  calculate_reaction_similarity(input_smiles, value)
# #         if result[similarity_type] >= cutoff_score:
# #             results_list.append(result)
# #     return results_list
# def get_top_similar_reactions_with_ids(file, input_smiles, cutoff_score, similarity_type):
#     values_list = parse_data(file)
#     results_list = []
    
#     # Perform similarity search for each value in values_list
#     for reaction_id, value in zip(get_reaction_ids(file), values_list):
#         result = calculate_reaction_similarity(input_smiles, value)
        
#         if result[similarity_type] >= cutoff_score:
#             results_list.append((reaction_id, result))
    
#     return results_list

# def get_top_similar_reactions_with_param(file, input_smiles, cutoff_score, similarity_type, params):
#     values_list = parse_data(file)
#     params_result_list = []
    
#     # Perform similarity search for each value in values_list
   
#     for parameter, value in zip(get_reaction_param(file,params), values_list):
#         result = calculate_reaction_similarity(input_smiles, value)
        
#         if result[similarity_type] >= cutoff_score:
#             params_result_list.append((parameter, result))
#             params_result_list = sorted(params_result_list, key=lambda x: x[1][similarity_type], reverse=True)
    
#     return params_result_list



# def get_rxn_info(param, reaction_id, file, input_smiles, cutoff_score, similarity_type):
#     param_list = param 
#     similar_reactions = get_top_similar_reactions_with_ids(file, input_smiles, cutoff_score, similarity_type)
#     for reaction in similar_reactions:
#         for reaction[0] in reaction:


input_smiles = '[CH3:1][NH:2][CH3:3].Cl[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1>>[CH3:1][N:2]([CH3:3])[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1'

file = '/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb'
similarity_type= 'total similarity'
cutoff_score = 0.5
params = ['identifiers', 'reaction_id']
reaction_similarity = ReactionSimilarity(file, input_smiles, cutoff_score, similarity_type, params)

result= reaction_similarity.get_top_similar_reactions_with_param()
print(result)
#  if mol1 is not None and mol2 is not None:
#         fp1 = AllChem.GetMorganFingerprint(mol1, 2)
#         fp2 = AllChem.GetMorganFingerprint(mol2, 2)
#         similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
#         return similarity
#     else:
#         return None

# # Now, 'similarities' contains the similarity scores for each value in values_list
# print(similarities)

    # dataset = dataset_pb2.Dataset()
    # dataset.ParseFromString(data)

    # result = {}

    # result['Dataset Name'] = dataset.name
    # result['Dataset Description'] = dataset.description

    # reactions_info = []
    # for reaction in dataset.reactions:
    #     reaction_identifiers = reaction.identifiers
    #     reaction_identifiers_list = []
    #     for identifier in reaction_identifiers:
    #         identifier_dict = {
    #             'type': identifier.type, 
    #             'value': identifier.value,
    #             'is_mapped': identifier.is_mapped
               
    #         }
    #         reaction_identifiers_list.append(identifier_dict)
    #     reactions_info.append(reaction_identifiers_list)

    # result['Reactions'] = reactions_info

    # return result



    


# result = parse_data(file_path)

# print(result)


# @app.route('/get_parse_data', methods=['POST'])
# def get_parse_data():
#     try:
#         data = request.get_json()
#         file_path = data.get('file_path')
#         rxn_data = parse_data(file_path)

#         if rxn_data is not None:
#             return jsonify(rxn_data), 200
#         else:
#             return jsonify({'error': 'Failed to get data'}), 400

#     except Exception as e:
#         print(f"Exception: {e}")
#         return jsonify({'error': str(e)}), 400

# if __name__ == "__main__":
#     app.run(port=7000)
     

# curl -X POST -H "Content-Type: application/json" -d '{"filepath": '/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb'}' http://localhost:7000/get_parse_data

# curl -X POST -H "Content-Type: application/json" -d "{\"file_path\": \"/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb\"}" http://localhost:7000/get_parse_data

# #     reaction_dict = {
# #         'reaction_id': reaction.reaction_id,
# #         'reaction_input': {
# #             'input_compound': {
# #                 'compound': {
# #                     'compound_id': reaction.reaction_input.input_compound.compound.compound_id,
# #                     'smiles': reaction.reaction_input.input_compound.compound.smiles,
# #                     'inchi': reaction.reaction_input.input_compound.compound.inchi,
# #                     # Add other compound fields as needed
# #                 },
# #                 'amount': {
# #                     'value': reaction.reaction_input.input_compound.amount.value,
# #                     'unit': reaction.reaction_input.input_compound.amount.unit,
# #                 },
# #                 'addition_device': {
# #                     'type': reaction.reaction_input.input_compound.addition_device.type,
# #                     'description': reaction.reaction_input.input_compound.addition_device.description,
# #                     # Add other addition_device fields as needed
# #                 },
# #                 # Add other input_compound fields as needed
# #             },
# #             'reaction_conditions': {
# #                 'temperature': {
# #                     'value': reaction.reaction_input.reaction_conditions.temperature.value,
# #                     'unit': reaction.reaction_input.reaction_conditions.temperature.unit,
# #                 },
# #                 'pressure': {
# #                     'value': reaction.reaction_input.reaction_conditions.pressure.value,
# #                     'unit': reaction.reaction_input.reaction_conditions.pressure.unit,
# #                 },
# #                 'stirring_conditions': {
# #                     'stirring_rate': {
# #                         'value': reaction.reaction_input.reaction_conditions.stirring_conditions.stirring_rate.value,
# #                         'unit': reaction.reaction_input.reaction_conditions.stirring_conditions.stirring_rate.unit,
# #                     },
# #                     'stirring_method': {
# #                         'type': reaction.reaction_input.reaction_conditions.stirring_conditions.stirring_method.type,
# #                         'description': reaction.reaction_input.reaction_conditions.stirring_conditions.stirring_method.description,
# #                         # Add other stirring_method fields as needed
# #                     },
# #                     # Add other stirring_conditions fields as needed
# #                 },
# #                 # Add other reaction_conditions fields as needed
# #             },
# #             # Add other input fields as needed
# #         },
# #         'reaction_outcome': {
# #             'products': [
# #                 {
# #                     'product_compound': {
# #                         'compound': {
# #                             'compound_id': product.product_compound.compound.compound_id,
# #                             'smiles': product.product_compound.compound.smiles,
# #                             'inchi': product.product_compound.compound.inchi,
# #                             # Add other compound fields as needed
# #                         },
# #                         'amount': {
# #                             'value': product.product_compound.amount.value,
# #                             'unit': product.product_compound.amount.unit,
# #                         },
# #                         'selectivity': {
# #                             'type': product.product_compound.selectivity.type,
# #                             'description': product.product_compound.selectivity.description,
# #                             # Add other selectivity fields as needed
# #                         },
# #                         # Add other product_compound fields as needed
# #                     },
# #                     'measurement': {
# #                         'mass_spec_measurement_details': {
# #                             'mass_spec_measurement_type': product.measurement.mass_spec_measurement_details.mass_spec_measurement_type,
# #                             # Add other mass_spec_measurement_details fields as needed
# #                         },
# #                         # Add other measurement fields as needed
# #                     },
# #                     # Add other product fields as needed
# #                 }
# #                 for product in reaction.reaction_outcome.products
# #             ],
# #             'analytical_data': [
# #                 {
# #                     'analyte': {
# #                         'compound': {
# #                             'compound_id': data.analyte.compound.compound_id,
# #                             'smiles': data.analyte.compound.smiles,
# #                             'inchi': data.analyte.compound.inchi,
# #                             # Add other compound fields as needed
# #                         },
# #                         'amount': {
# #                             'value': data.analyte.amount.value,
# #                             'unit': data.analyte.amount.unit,
# #                         },
# #                         # Add other analyte fields as needed
# #                     },
# #                     # Add other analytical_data fields as needed
# #                 }
# #                 for data in reaction.reaction_outcome.analytical_data
# #             ],
# #             # Add other reaction_outcome fields as needed
# #         },
# #         # Add other top-level fields as needed
# #     }
