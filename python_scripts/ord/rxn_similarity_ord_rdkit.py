# import dataset_pb2
# from rdkit import Chem
# from rdkit.Chem import AllChem

# def parse_data(file_path):
#     with open(file_path, 'rb') as f:
#         data = f.read()

#     dataset = dataset_pb2.Dataset()
#     dataset.ParseFromString(data)

#     result = {}

#     result['Dataset Name'] = dataset.name
#     result['Dataset Description'] = dataset.description

#     reactions_info = []
#     input_smarts = '...'  # Replace with your input SMARTS string

#     # Convert the input SMARTS string to an RDKit reaction object
#     input_reaction = AllChem.ReactionFromSmarts(input_smarts)

#     # Perform similarity search
#     similarity_threshold = 0.8  # Adjust as needed
#     for reaction in dataset.reactions:
#         # Convert the Protocol Buffers reaction to an RDKit reaction object
#         pb_reaction_smarts = reaction.identifiers
#         pb_reaction = AllChem.ReactionFromSmarts(pb_reaction_smarts)

#         # Perform similarity comparison
#         similarity = AllChem.GetBestRMS(input_reaction, pb_reaction)

#         # Include reactions that meet the similarity threshold
#         if similarity >= similarity_threshold:
#             reaction_info = {
#                 'Reaction ID': reaction.reaction_id,
#                 'Similarity': similarity,
#                 'Reaction Identifier': pb_reaction_smarts
#             }
#             reactions_info.append(reaction_info)

#     result['Reactions'] = reactions_info

#     return result

# file_path = '/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb'
# result = parse_data(file_path)

# print(result)

from ord_schema.proto import dataset_pb2, reaction_pb2
from google.protobuf import message
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_data(file_path):
    with open(file_path, 'rb') as f:
        data = f.read()

    dataset = dataset_pb2.Dataset()
    dataset.ParseFromString(data)

    result = {}

    result['Dataset Name'] = dataset.name
    result['Dataset Description'] = dataset.description

    reactions_info = []
    input_smarts = '[Cl:1][C:2]1[CH:12]=[N:11][C:10]2[O:9][CH:8]([CH2:13][CH2:14]Cl)[CH2:7][N:6]([CH3:16])[C:5](=[S:17])[C:4]=2[CH:3]=1.[C:18](=[O:21])([O-:20])[O-].[K+].[K+].[NH:24]1[CH2:27][CH2:26][CH2:25]1.CS(C)=[O:30]>O>[C:10]([OH:9])(=[O:30])[C:18]([OH:20])=[O:21].[N:24]1([CH2:14][CH2:13][CH:8]2[CH2:7][N:6]([CH3:16])[C:5](=[S:17])[C:4]3[CH:3]=[C:2]([Cl:1])[CH:12]=[N:11][C:10]=3[O:9]2)[CH2:27][CH2:26][CH2:25]1 |f:1.2.3,7.8|'

    # Convert the input SMARTS string to an RDKit reaction object
    input_reaction = AllChem.ReactionFromSmarts(input_smarts)

    # Perform similarity search
    similarity_threshold = 0.8  # Adjust as needed
    for reaction in dataset.reactions:
        # Find the reaction input SMILES in the repeated field
        smarts_message = next(smiles for smiles in reaction.identifiers if smiles.type == dataset_pb2.Reaction.ReactionIdentifier.REACTION_CXSMILES)
        
        # Extract the value (SMILES/SMARTS) from the message
        pb_reaction_smarts = smarts_message.value

        # Convert the Protocol Buffers reaction to an RDKit reaction object
        pb_reaction = AllChem.ReactionFromSmarts(pb_reaction_smarts)

        # Perform similarity comparison
        similarity = AllChem.GetBestRMS(input_reaction, pb_reaction)

        # Include reactions that meet the similarity threshold
        if similarity >= similarity_threshold:
            reaction_info = {
                'Reaction ID': reaction.reaction_id,
                'Similarity': similarity,
                'Reaction Input': pb_reaction_smarts
            }
            reactions_info.append(reaction_info)

    result['Reactions'] = reactions_info

    return result

file_path = '/Users/ashort/ord_data/ord-data/data/01/ord_dataset-01dbb772c5e249108f0b191ed17a2c0c.pb'
result = parse_data(file_path)

print(result)



