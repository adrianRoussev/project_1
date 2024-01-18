from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS


def get_mols_from_rxn(rxn_smarts):

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

# def index_reactants_or_products(mol):
#     total_atom_count = 0
#     for atom in mol.GetAtoms():
#         atom.SetAtomMapNum(atom.GetIdx() + 1 + total_atom_count) 
#         total_atom_count += mol.GetNumAtoms() 



def get_rxn_from_mols(rxn_dict):
    reactants = [Chem.MolFromSmiles(reactant) for reactant in rxn_dict['reactants']]
    products = [Chem.MolFromSmiles(product) for product in rxn_dict['products']]
    catalysts = [Chem.MolFromSmiles(catalyst) for catalyst in rxn_dict['catalysts']]

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

    
    rxn_smarts = AllChem.ReactionToSmarts(reaction)
    return rxn_smarts

def get_indexed_smarts(smiles_string):

    mol = Chem.MolFromSmiles(smiles_string)
   
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx()+1)
    mol

    smarts_with_indices = Chem.MolToSmarts(mol)

    print(smarts_with_indices)




def highlight_similarity(rxn_dict):
    reactants = [Chem.MolFromSmiles(reactant) for reactant in rxn_dict['reactants']]
    products = [Chem.MolFromSmiles(product) for product in rxn_dict['products']]
   
    

    results_dict = {}
    smarts_product_pair = {}
    product_atom_pair = {}
    reactant_atoms_with_match = []
    index = 0
    atoms_with_match_mol = {}
    for reactant in reactants:
        smarts_list = []
       
        for product in products:
            mols = [reactant,product]
            result =rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny)
            result_smarts = result.smartsString
            smarts_product_pair[result_smarts] = product
            smarts_list.append(result_smarts)
        smarts_list_sorted = sorted(smarts_list, key=len, reverse=True)
        match = smarts_list_sorted[0]
        pattern = Chem.MolFromSmarts(match)
        
        reactant_atoms = list(reactant.GetSubstructMatch(pattern))
        
        match_product = smarts_product_pair[match]
        product_atoms = list(match_product.GetSubstructMatch(pattern))
        
        match_product_smarts = Chem.MolToSmarts(match_product)
  
       
        for atom in reactant_atoms:
            atom_index = atom + index + 1 
            if atom_index not in atoms_with_match_mol:
                reactant_atoms_with_match.append(atom_index)
        index += reactant.GetNumAtoms()

        results_dict[Chem.MolToSmarts(reactant)]= {'reactant atoms': reactant_atoms, 
                                                    'match product': match_product_smarts, 
                                                    'matching smarts': match, 
                                                    'product atoms':product_atoms}
    reactant_atom_count = 0
    product_atom_count = 0 
    
    for reactant_mol in reactants: 
        
        reactant_mol_atom_count = reactant_mol.GetNumAtoms()
        reactant_atom_count += reactant_mol_atom_count

    for  product_mol in  products: 

        product_mol_atom_count = product_mol.GetNumAtoms()
        product_atom_count += product_mol_atom_count

    match_reactant_atom_count = 0

    for atom in reactant_atoms_with_match:
        match_reactant_atom_count += 1

    
    expected_atoms = list(range(1,reactant_atom_count))

    atoms_without_match_idx = []
    atoms_without_match = []
    if reactant_atoms_with_match != expected_atoms:
        for atom in expected_atoms:
            if atom not in reactant_atoms_with_match: 
                atoms_without_match_idx.append(atom)

    
    missing_patterns = {}
    missing_product_patterns = {}
    product_atoms_with_match = []
    product_atoms_without_match = []
    if atoms_without_match_idx:
        numbers = []
    
        number_of_atoms = 0
        for reactant in reactants:
            atom_group = {}
            idxs = []
            atomic_nums = []
            for atom in atoms_without_match_idx:
                number= atom-1-number_of_atoms
                
                if  number_of_atoms < atom <= reactant.GetNumAtoms() + number_of_atoms:
                    atom_with_idx = reactant.GetAtomWithIdx(number)
        
                    atomic_num = atom_with_idx.GetAtomicNum()
                    atoms_without_match.append([atomic_num,atom])
                    idxs.append(atom)
                    atomic_nums.append(atomic_num)
            if idxs:
                atom_group['indecies']= idxs
                atom_group['atomic numbers'] = atomic_nums
            if atom_group:
                missing_patterns[Chem.MolToSmarts(reactant)] = atom_group
                
                matching_product_atoms = results_dict[Chem.MolToSmarts(reactant)]['product atoms']
                matching_product_atoms_count = 0 
                for each in matching_product_atoms:
                   matching_product_atoms_count += 1

                product_atom_index = matching_product_atoms_count + number_of_atoms + 1
                matching_product = results_dict[Chem.MolToSmarts(reactant)]['match product']
                # matching_product_index = products.index(Chem.MolFromSmarts(matching_product))
                product_atoms_with_match.append(product_atom_index)
                
            number_of_atoms += reactant.GetNumAtoms()

        expected_product_atoms = list(range(1,product_atom_count))
        for atoms in expected_product_atoms:    
            if atoms not in product_atoms_with_match:
            
                product_atoms_without_match.append(atoms)

        product_atom_increments = []
        
        product_atom_count = 0

        for product in products: 
            count =  product.GetNumAtoms()
            product_atom_count += count
            product_atom_increments.append(product_atom_count)
       
        
        for increment in product_atom_increments: 
            
            atom_groups_product = {}
            patterns_product = []
            idxs_product = []
           

            index = product_atom_increments.index(increment)
            product = products[index]
            previous_index = index - 1 
            previous_increment = product_atom_increments[previous_index]

            pattern = []
            pattern_idx = []

            for atom in product_atoms_without_match: 
                 
                if 0 <= previous_index  :
                    previous_increment = product_atom_increments[previous_index]
                else:
                    previous_increment= 0  
                product_index = atom -1 - previous_increment
                atomic_num = product.GetAtomWithIdx(product_index) 
                    
                if previous_increment < atom <= increment:

                    if product_index == 0 or pattern_idx[-1]  == product_index - 1:
                        pattern.append(Chem.SmartsFromAtom(atomic_num))
                        pattern_idx.append(product_index)
                        
                    if atom == increment:
                        idxs_product.append(pattern_idx)
                        patterns_product.append(pattern)
                    
                    if pattern_idx[-1]  != product_index - 1:
                        idxs_product.append(pattern_idx)
                        patterns_product.append(pattern)
                        pattern_idx = [product_index]
                        pattern = [Chem.SmartsFromAtom(atomic_num)]

            atom_groups_product['product'] = Chem.MolToSmarts(product)
            atom_groups_product['patterns'] = patterns_product
            atom_groups_product['idxs'] = idxs_product



                        

        # for increment in product_atom_increments:
        #     for atom in product_atoms_without_match:
                
        #         if product_atom_increments[product_atom_increments.index(increment) - 1] < atom <= increment:
        #             product_index = product_atom_increments.index(increment)






    print(results_dict,reactant_atom_count,  product_atom_count, match_reactant_atom_count, atoms_without_match, atoms_without_match_idx,  reactant_atoms_with_match, missing_patterns, product_atoms_without_match, product_atom_increments, atom_groups_product) 
#IS SUBSTRUCTURE MATCH NOT COUNTING ATOM MATCH IF BOND IS different



smiles_string = "CCO"
reactants = ['CCO', 'CC(=O)O']
products = ['CC(=O)OC', 'CCO']
agents = ['[H+]', '[Na+]']
rxn_smarts= "[Mg:1].[C:2](=[O:3])>>[C:2]([Mg:1])[OH:3]"

# 1. Aldol Condensation:
#    [C:1](=[O:2])-[CH2:3]-[C:4](=[O:5])>>[C:1](=[O:2])-[C:3]([OH:6])-[C:4](=[O:5])

# 2. Grignard Reaction:
#    [Mg:1].[C:2](=[O:3])>>[C:2]([Mg:1])[OH2:3]
#[Mg][CX3:1][X:2].[C:3](=[O:4])[H:5]>>[C:3][X:2][CX3:1][O:4][H:5]
# 3. Ring Closure - Diels-Alder:
#    [C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3]1=[C:4]

# 4. Substitution - SN2:
#    [C:1]-[Br:2].[N:3]>>[C:1]-[N:3].[Br:2]

# 5. Oxidation of Alcohol:
#    [C:1]-[OH:2]>>[C:1](=[O:3])-[OH:2]

# 6. Reduction of Ketone:
#    [C:1](=[O:2])-[H:3]>>[C:1](-[H:3])-[OH1-:2]

# 'CCO.CC(=O)O>>CC(=O)OC.CCO'

rxn_dict = get_mols_from_rxn(rxn_smarts)
smarts= get_rxn_from_mols(rxn_dict)
mol = get_indexed_smarts(smiles_string)
highlight_similarity(rxn_dict)
print("Reaction SMARTS:", smarts)

# # Example usage:
# rxn_smarts_example = 
# rxn_dict_example = get_mols_from_rxn(rxn_smarts_example)
# print("Reaction Dictionary:", rxn_dict_example)

# rxn_smarts_reconstructed = get_rxn_from_mols(rxn_dict_example)
# print("Reconstructed Reaction SMARTS:", rxn_smarts_reconstructed)


#  def get_rxn_from_mols(rxn_dict):
  
#     reactants = [Chem.MolFromSmiles(reactant) for reactant in rxn_dict['reactants']]
#     products = [Chem.MolFromSmiles(product) for product in rxn_dict['products']]
#     catalysts = [Chem.MolFromSmiles(catalyst) for catalyst in rxn_dict['catalysts']]

#     molecule_list = []

#     for mol in reactants:
#         mol.SetProp("_RXN_ROLE", "reactant")
#         molecule_list.append(mol)


#     for mol in catalysts:
#         mol.SetProp("_RXN_ROLE", "agent")
#         molecule_list.append(mol)

    
#     for mol in products:
#         mol.SetProp("_RXN_ROLE", "product")
#         molecule_list.append(mol)

#     molecules = (molecule_list)

#     reaction = AllChem.ReactionFromMolecule(molecules)
    
#     rxn_smarts = AllChem.ReactionToSmarts(reaction)
#     return rxn_smarts