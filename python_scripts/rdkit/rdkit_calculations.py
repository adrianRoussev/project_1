from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS
from itertools import permutations
from itertools import product
from itertools import combinations


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


# def find_closest_match(reactant, product_list):
    
#     best_match = None
#     best_mcs = None

#     for bond_compare in [rdFMCS.BondCompare.CompareOrderExact, rdFMCS.BondCompare.CompareOrder, rdFMCS.BondCompare.CompareAny]:
#         for atom_compare in [rdFMCS.AtomCompare.CompareIsotopes, rdFMCS.AtomCompare.CompareElements, rdFMCS.AtomCompare.CompareAny]:
#             for complete_rings_only in [True, False]:
#                 mcs_params = rdFMCS.MCSParameters()
#                 mcs_params.BondCompare = bond_compare
#                 mcs_params.AtomCompare = atom_compare
#                 mcs_params.CompleteRingsOnly = complete_rings_only

#                 for product in product_list:
                   

#                     # Find MCS
#                     mcs = rdFMCS.FindMCS([reactant, product], mcs_params)

#                     if mcs.numAtoms and (best_mcs is None or mcs.numAtoms > best_mcs.numAtoms):
#                         best_mcs = mcs
#                         best_match = product
#     results = {}

#     results['best_msc'] = best_mcs.smartsString
#     results['product_match'] = best_match

#     return results

# def find_closest_match(reactant_smarts, product_list):

#     best_matches = []
    # best_mcs_list = []

    # for bond_compare in [rdFMCS.BondCompare.CompareOrderExact, rdFMCS.BondCompare.CompareOrder, rdFMCS.BondCompare.CompareAny]:
    #     for atom_compare in [rdFMCS.AtomCompare.CompareIsotopes, rdFMCS.AtomCompare.CompareElements, rdFMCS.AtomCompare.CompareAny]:
    #         for complete_rings_only in [True, False]:
    #             mcs_params = rdFMCS.MCSParameters()
    #             mcs_params.BondCompare = bond_compare
    #             mcs_params.AtomCompare = atom_compare
    #             mcs_params.CompleteRingsOnly = complete_rings_only

    #             matches = []

    #             for product in product_list:
    #                 # Find MCS

    #                 mcs = rdFMCS.FindMCS([reactant, product], mcs_params)

    #                 if mcs.numAtoms:
    #                     matches.append((mcs.numAtoms, product, mcs))

    #             if matches:
    #                 matches.sort(reverse=True, key=lambda x: x[0])
    #                 best_matches.append(matches[0][1])
    #                 best_mcs_list.append(matches[0][2])

    # results = {}

    # results['best_mcs_list'] = [mcs.smartsString for mcs in best_mcs_list]
    # results['product_matches'] = best_matches

    # return results

def indexes_of_matching_elements(main_list, sub_list):
    indexes = []
    sub_list_length = len(sub_list)
    for i in range(len(main_list) - sub_list_length + 1):
        if main_list[i:i+sub_list_length] == sub_list:
            indexes.append(i)
     
    return indexes


def find_closest_match(reactant, product_list):
    best_matches = []
    best_mcs_list = []
    matches = {}

    for product in product_list:

        params = rdFMCS.MCSParameters()
        params.AtomCompare = rdFMCS.AtomCompare.CompareElements
        params.BondCompare = rdFMCS.BondCompare.CompareAny
        params.RingMatchesRingOnly = False

        all_atoms_idxs = (set(range(product.GetNumAtoms())))
           
        # Find MCS for the first check
        mcs1 = rdFMCS.FindMCS([reactant, product], params)
        matched_atoms = Chem.MolFromSmarts(mcs1.smartsString)

                
        
        matched_atoms_idxs_list = product.GetSubstructMatches(matched_atoms)
        matched_atoms_idxs = []
        
        for each in matched_atoms_idxs_list:
            for atom in each:
                matched_atoms_idxs.append(atom)
        # Retain the first and last matched atoms, and all unmatched atoms
        atoms_to_retain = set((all_atoms_idxs - set(matched_atoms_idxs)))

        # Create a new molecule and remove atoms not in atoms_to_retain
        product_retained_atoms = Chem.RWMol(product)
        for atom_idx in reversed(range(product.GetNumAtoms())):
            if atom_idx not in atoms_to_retain:
                product_retained_atoms.RemoveAtom(atom_idx)

        product_retained_atoms = product_retained_atoms.GetMol()

        # Find MCS for the second check with the remaining product
        remaining_atoms= product_retained_atoms.GetAtoms()
        reactant_atoms = reactant.GetAtoms()

        remaining_atom_list = []
        reactant_atom_list = []
        
        for each in remaining_atoms:
            remaining_atom_list.append(each.GetAtomicNum())
                    
        for each in reactant_atoms:
            reactant_atom_list.append(each.GetAtomicNum())
        
        matching_substruct = indexes_of_matching_elements(remaining_atom_list,reactant_atom_list)
        index_sets = []
        for each in matching_substruct:
            index_set = []
            count = 0
            for num in reactant_atom_list:
                index_set.append(each+count)
                count += 1
            index_sets.append(index_set)

        total_matched_atoms = []

        for atoms in matched_atoms_idxs_list:
                group = []
                for atom in atoms:
                    group.append(atom)
                total_matched_atoms.append(group)
        
        for index in index_sets: 
            total_matched_atoms.append(index)

        
        matches[Chem.MolToSmarts(product)] = total_matched_atoms
                    
    return matches


                    # if mcs.numAtoms:
                        # matches.append((mcs.numAtoms, Chem.MolToSmarts(product), mcs.smartsString))

    #                     best_matches.append(matches[0][1])
    #             #     best_mcs_list.append(matches[0][2])
    
   
    # print(matches)

    
    # return best_mcs_list, best_matches


# def find_non_overlapping_matches(reactants, products):
#     result_dict = {}

#     for reactant in reactants:
#         reactant_index = reactants.index(reactant)
#         matches_info = []

#         # Use the find_closest_match function to get initial matches
#         best_mcs_list, best_matches = find_closest_match(reactant, products)

#         for best_mcs, best_match in zip(best_mcs_list, best_matches):
#             smarts= best_mcs.smartsString
#             mol =Chem.MolFromSmarts(smarts)
#             reactant_match_atoms = set(reactant.GetSubstructMatch(mol))
#             product_match_atoms = set(reactant_match_atoms).union({len(reactant.GetAtoms())})

#             overlap = False
#             for existing_match_info in matches_info:
#                 if reactant_match_atoms & existing_match_info['reactant_match_atoms']:
#                     overlap = True
#                     break

#     if not overlap:
#         match_info = {
#             'product': best_match,
#             'matching_substructure': best_mcs.smartsString,
#             'reactant_match_atoms': reactant_match_atoms,
#             'product_match_atoms': product_match_atoms,
#         #     'reactant_atom_indices': best_mcs.MCSAtomIndices[0],
#         #     'product_atom_indices': [i + len(reactant.GetAtoms()) for i in best_mcs.MCSAtomIndices[1]],
#         }

#         matches_info.append(match_info)

#     result_dict[reactant] = matches_info

                
#     return result_dict



def find_non_overlapping_matches(reactants, products):
    result_dict = {}

    for reactant in reactants:
        reactant_index = reactants.index(reactant)
        matches_info = []

        # Use the find_closest_match function to get initial matches
        best_mcs_list, best_matches = find_closest_match(reactant, products)

        for best_mcs, best_match in zip(best_mcs_list, best_matches):
            smarts = best_mcs.smartsString
            mol = Chem.MolFromSmarts(smarts)
            reactant_match_atoms = set(reactant.GetSubstructMatch(mol))
            product_match_atoms = set(reactant_match_atoms).union({len(reactant.GetAtoms())})

            overlap_count = 0
            for existing_match_info in matches_info:
                common_atoms = reactant_match_atoms & existing_match_info['reactant_match_atoms']
                if len(common_atoms) > 1:
                    overlap_count += len(common_atoms)
                    break

            if overlap_count <= 1:
                match_info = {
                    'product': Chem.MolToSmarts(best_match),
                    'matching_substructure': best_mcs.smartsString,
                    'reactant_match_atoms': reactant_match_atoms,
                    'product_match_atoms': product_match_atoms,
                    # 'reactant_atom_indices': best_mcs.MCSAtomIndices[0],
                    # 'product_atom_indices': [i + len(reactant.GetAtoms()) for i in best_mcs.MCSAtomIndices[1]],
                }

                matches_info.append(match_info)

        result_dict[Chem.MolToSmarts(reactant)] = matches_info

    return result_dict

# def find_best_combination(matches_dict):
#     # Get the list of reactants and their matches
#     reactants = list(matches_dict.keys())
#     matches_info_list = list(matches_dict.values())

#     # Generate all possible permutations of matches
#     all_permutations = list(permutations(matches_info_list))

#     best_combination = None
#     max_matched_atoms = 0

#     # Iterate through each permutation
#     for permutation in all_permutations:
#         combined_matches = []
#         matched_atoms = set()

#         # Combine matches from different reactants
#         for reactant_matches in permutation:
#             for match_info in reactant_matches:
#                 product_atoms = match_info['product_match_atoms']
#                 if not (matched_atoms & product_atoms):
#                     combined_matches.append(match_info)
#                     matched_atoms |= product_atoms

#         # Check if the current combination covers more atoms than the previous best
#         if len(matched_atoms) > max_matched_atoms:
#             max_matched_atoms = len(matched_atoms)
#             best_combination = combined_matches

#     return best_combination

# def find_best_combination_with_overlap(reactants, products):
#     matches_dict = {}

#     # Obtain non-overlapping matches for each reactant
#     for reactant in reactants:
#         matches_dict[reactant] = find_non_overlapping_matches([reactant], products)[Chem.MolToSmarts(reactant)]

#     # Find the best combination of matches
#     best_combination = find_best_combination(matches_dict)

#     return best_combination


# def find_best_combination(matches_dict):
#     # Get the list of reactants and their matches
#     reactants = list(matches_dict.keys())
#     matches_info_lists = list(matches_dict.values())

#     # Generate all possible combinations of matches from different reactants
#     all_combinations = list(product(*matches_info_lists))

#     best_combination = None
#     max_matched_atoms = 0

#     # Iterate through each combination
#     for combination in all_combinations:
#         combined_matches = []
#         matched_atoms = set()

#         # Combine matches from different reactants
#         for reactant_matches in combination:
#             for match_info in reactant_matches:
#                 product_atoms = match_info['product_match_atoms']
#                 if not (matched_atoms & product_atoms):
#                     combined_matches.append(match_info)
#                     matched_atoms |= product_atoms

#         # Check if the current combination covers more atoms than the previous best
#         if len(matched_atoms) > max_matched_atoms:
#             max_matched_atoms = len(matched_atoms)
#             best_combination = combined_matches

#     return best_combination


# def find_best_combination(reactants, products):
#     matches_dict = {}

#     # Get the matches_dict for each reactant
#     for reactant in reactants:
#         matches_dict[reactant] = find_closest_match(reactant, products)

#     # Get the list of reactants and their matches
#     reactants = list(matches_dict.keys())
#     matches_info_lists = list(matches_dict.values())

#     # Generate all possible combinations of matches from different reactants
#     all_combinations = list(product(*matches_info_lists))

#     best_combination = None
#     max_matched_atoms = 0

#     # Iterate through each combination
#     for combination in all_combinations:
#         combined_matches = []
#         matched_atoms = set()

#         # Combine matches from different reactants
#         for reactant_matches in combination:
#             for match_info in reactant_matches:
#                 product_atoms_result = (match_info['product_match_atoms'])
#                 product_atoms = product_atoms_result.smartsString
#                 if not (matched_atoms & product_atoms):
#                     combined_matches.append(match_info)
#                     matched_atoms |= product_atoms

#         # Check if the current combination covers more atoms than the previous best
#         if len(matched_atoms) > max_matched_atoms:
#             max_matched_atoms = len(matched_atoms)
#             best_combination = combined_matches

#     return best_combination

def find_value_combination(input_set, input_list):
   
    total_combinations = []
    for item in input_list:
        for pattern, combos in item.items():
            for combination in combos:
                total_combinations.append(combination)
    # print(total_combinations)
    count = 0
    for each in total_combinations:
        count+=1
    
    index_set = set(range(1,count))
    
    unique_pairs = combinations(index_set,2)
    
    print(unique_pairs)
    
    list_combinations = []
    for pair in unique_pairs:
        list_pair = []
        for index in pair:
            list_pair.append(total_combinations[index])
        list_combinations.append(list_pair)

    print(list_combinations)
    
    value_list = []
    for combination in list_combinations:
        value_set = []
        for list in combination:
            for value in list:
                value_set.append(value)
        value_list.append(value_set)
    
    missing_values_list = []
    for value_set in value_list:
        missing_values = []
        for value in input_set:
            if value not in value_set:
                missing_values.append(value)
        missing_values_list.append(missing_values)

    print(missing_values_list, ":)")
    
    no_missing_values_index = []
    missing_values_dict ={}
    
    for each in missing_values_list:
        if each == []:
            count = 0
            no_missing_values_index.append(missing_values_list.index(each))
    
        if no_missing_values_index == []:
        # for each in missing_values_list:
            count = each.index(each[-1]) + 1
            missing_values_dict[count] = [each, missing_values_list.index(each)]
            # missing_values_count.append = count
            # if missing_values == []:
            #     missing_values.append = each
            # elif missing_values[-1] == count
    print(missing_values_dict)
    keys_list = []
    for each in missing_values_dict.keys():
        keys_list.append(each)
    if keys_list: #remove after 'if no_missing_values_index == []:' indenting is fixed
        min_key = min(keys_list)

        min_dict = {key: value for key, value in missing_values_dict.items() if key == min_key}

        print(min_dict)

        
        
    # sorted_dict = dict(sorted( missing_values_dict.items(), key=lambda item: item[1]))
    # print(missing_values_dict)
    # selected_pairs= {}
    # lowest_value = []
    # for key,value in sorted_dict:
    #     if selected_pairs == {}:
    #         selected_pairs[key] = value 
    #         lowest_value.append(value)
    #     else:
    #         if value == lowest_value:
    #             selected_pairs[key] = value
    # print(selected_pairs)
    # selected = []
    # for key, value in selected_pairs:
    #     if value in key == input_set[0] or input_set[-1]:

    selected_combinations = []
    if no_missing_values_index:
        for index in no_missing_values_index:
            selected_combinations.append(list_combinations[index])
    elif min_dict:
        if len(min_dict) ==1:
            selected_combinations.append(min_dict.values()[0])
        else:
            #combinations that are favourable 
    
    print(selected_combinations)

    # min_overlap = float('inf')
    # selected_elements = []

    # for elements in total_combinations:
    #     covered_values = set(elements)
    #     overlap_count = len(covered_values & input_set)

    #     if overlap_count < min_overlap and covered_values == input_set:
    #         min_overlap = overlap_count
    #         selected_elements = [elements]
    #     elif overlap_count == min_overlap and covered_values == input_set:
    #         selected_elements.append(elements)
    # print(selected_elements)
    # return selected_elements
    # permuted_combinations = permutations(total_combinations)
    # # unique_permutations = set(permuted_combinations)
    # print(permuted_combinations)
    # # count_unique_elements = sum(1 for perm in unique_permutations if set(perm) == input_set)
    # print(count_unique_elements)
               
    # if count_unique_elements > best_count:
    #     best_combination = combination
    #     best_count = count_unique_elements
    # print(best_combination, best_count)
    # return best_combination


def select_elements(group, target_set):
    
    min_overlap = float('inf')
    selected_elements = []

    for elements in group:
        covered_values = set(elements)
        overlap_count = len(covered_values & target_set)

        if overlap_count < min_overlap and covered_values == target_set:
            min_overlap = overlap_count
            selected_elements = [elements]
        elif overlap_count == min_overlap and covered_values == target_set:
            selected_elements.append(elements)

    return selected_elements
def find_best_combination(reactants, products):
    matches_dict = {}
    product_combos= {}
    product_result = {}
    for reactant in reactants:
        
        matches_dict[Chem.MolToSmarts(reactant)] = find_closest_match(reactant, products)

    for product in products:
        product_combos_list = []
        count = 0
        product_atom_idxs = []
        for atom in product.GetAtoms():
            product_atom_idxs.append(count)
            count+=1
        # product_atom_idxs = list(range(product.GetAtoms))
       
        print(product_atom_idxs)
        for reactant in reactants:
            combos = matches_dict[Chem.MolToSmarts(reactant)][Chem.MolToSmarts(product)]
            
            product_combos_list.append({Chem.MolToSmarts(reactant): matches_dict[Chem.MolToSmarts(reactant)][Chem.MolToSmarts(product)]})
        product_combos[Chem.MolToSmarts(product)] = product_combos_list
        result = find_value_combination(product_atom_idxs, product_combos_list)
        product_result[Chem.MolToSmarts(product)] = result

    print(product_result)


    # fix inputed range in above function so that you can continue here  


# def find_best_combination(reactants, products):
#     matches_dict = {}
#     total_match_combinations =  {}
#     # Get the matches_dict for each reactant
#     for reactant in reactants:
#         best_mcs_list, best_matches = find_closest_match(reactant, products)
#         matches_info = []

#         for best_mcs, best_match in zip(best_mcs_list, best_matches):
#             smarts = best_mcs.smartsString
#             mol = Chem.MolFromSmarts(smarts)
#             reactant_match_atoms = set(reactant.GetSubstructMatches(mol))
#             product_match_atoms = set(best_match.GetSubstructMatches(mol))
#             print(Chem.MolToSmarts(best_match))
#             # product_match_atoms = set(reactant_match_atoms).union({len(reactant.GetAtoms())})

#             overlap_count = 0
#             for existing_match_info in matches_info:
#                 common_atoms = reactant_match_atoms & existing_match_info['reactant_match_atoms']
#                 if len(common_atoms) > 1:
#                     overlap_count += len(common_atoms)
#                     break

#             # if overlap_count <= 1:
            
#             match_combos_reactant = []
#             match_info = {
#                     'product': Chem.MolToSmarts(best_match),
#                     'matching_substructure': best_mcs.smartsString,
#                     'reactant_match_atoms': reactant_match_atoms,
#                     'product_match_atoms': product_match_atoms,
#                 }
#                 # if [match] not in total_match_combinations:
#             total_match_combinations[reactant]= [product_match_atoms]
#             matches_info.append(match_info)
#             if reactant_match_atoms not in match_combos_reactant:
#                 match_combos_reactant.append(reactant_match_atoms)
#         matches_dict[Chem.MolToSmarts(reactant)] = matches_info
#     #from total_match_combinations for each reactant for each atom match make a list of all possible combinations between each group with at least one atom set from each group. find combo with closest total amount of atoms to total product atom count (equal or less than but not more atoms). also prioritize least overlap. 
#         # total_match_combinations[reactant]= match_combos_reactant

#     # all_combinations = list(product(*matches_dict.values()))
#     # # Get the list of reactants and their matches
#     # reactants = list(matches_dict.keys())
#     # matches_info_lists = list(matches_dict.values())
#     # matches_info_lists= {}
#     # for reactant in reactants:
#     #     matches_info_list= matches_dict[Chem.MolToSmarts(reactant)]
#     #     matches_info_lists[Chem.MolToSmarts(reactant)]= matches_info_list
#     # # Generate all possible combinations of matches from different reactants
#     all_combinations = list((total_match_combinations.values()))
#     # # for each reactant, if the same match product, then product_match_atoms can only have one in common 

#     # best_combination = None
#     # max_matched_atoms = 0

#     # # Iterate through each combination
#     # for combination in all_combinations:
#     #     combined_matches = []
#     #     matched_atoms = set()

#     #     # Combine matches from different reactants
#     #     for reactant_matches in combination:
#     #         for match_info in matches_info:
#     #             product_atoms = match_info['product_match_atoms']
#     #             reactant_atoms= match_info['reactant_match_atoms']
#     #             if not (matched_atoms & product_atoms):
#     #                 combined_matches.append(match_info)
#     #                 matched_atoms |= product_atoms

#     #     # # Check if the current combination covers more atoms than the previous best
#     #     # if len(matched_atoms) > max_matched_atoms:
#     #     #     max_matched_atoms = len(matched_atoms)
#     #     #     best_combination = combined_matches
#     print(all_combinations)
#     # return best_combination
    




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
    
        # for product in products:
            # mols = [reactant,product]
            # result =rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny)
        #     result_smarts = result.smartsString
        #     smarts_product_pair[result_smarts] = product
        #     smarts_list.append(result_smarts)
        # smarts_list_sorted = sorted(smarts_list, key=len, reverse=True)
        # match = smarts_list_sorted[0]
        result = find_closest_match(reactant, products)
        match = result['best_mcs_list']
        pattern = Chem.MolFromSmarts(match)
        
        reactant_atoms = list(reactant.GetSubstructMatch(pattern))
        
        # match_product = smarts_product_pair[match]
        match_product = result['product_matches']
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

    
    expected_atoms = list(range(1,reactant_atom_count+1))

    atoms_without_match_idx = []
    atoms_without_match = []
    if reactant_atoms_with_match != expected_atoms:
        for atom in expected_atoms:
            if atom not in reactant_atoms_with_match: 
                atoms_without_match_idx.append(atom)

    
    missing_patterns = {}
    missing_product_patterns = {}
    product_atoms_with_match = []
    product_atoms_with_match_atom = []
    product_atoms_without_match = []
    # if atoms_without_match_idx:
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
        if idxs != []:
            atom_group['indecies']= idxs
            atom_group['atomic numbers'] = atomic_nums
            # if atom_group != {}:
        missing_patterns[Chem.MolToSmarts(reactant)] = atom_group
                
        matching_product_atoms = results_dict[Chem.MolToSmarts(reactant)]['product atoms']
        matching_product_atoms_count = 0 
            #issue here... 
        for each in matching_product_atoms:
            matching_product_atoms_count += 1

            product_atom_index = each + number_of_atoms + 1
            matching_product = results_dict[Chem.MolToSmarts(reactant)]['match product']
                # matching_product_index = products.index(Chem.MolFromSmarts(matching_product))
            product_atoms_with_match.append(product_atom_index)
                #this has right index
            product_atoms_with_match_atom.append(each)
                
    number_of_atoms += reactant.GetNumAtoms()
    test = product_atom_count
    expected_product_atoms = list(range(1,test+1))
    for atoms in expected_product_atoms:    
        if atoms not in product_atoms_with_match:
            
            product_atoms_without_match.append(atoms)

    product_atom_increments = []
        
    product_atom_count = 0

    for product in products: 
        count =  product.GetNumAtoms()
        product_atom_count += count
        product_atom_increments.append(product_atom_count)
       
        # product_atom_increments = [3]
        # product_atoms_with_match = [2]
        # the overall rxn is "[Mg:1].[C:2](=[O:3])>>[C:2]([Mg:1])[OH:3]"
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
            
            if pattern_idx != []:
                last = pattern_idx[-1]
            else:
                last = "x"
                    
            if previous_increment < atom <= increment:
                
                product_index = atom -1 - previous_increment
                atomic_num = product.GetAtomWithIdx(product_index) 

                if last  != product_index - 1 and pattern_idx != []:
                
                    
                    idxs_product.append(pattern_idx)
                    patterns_product.append(pattern)
                    pattern_idx = [product_index]
                    pattern = [atomic_num.GetSmarts()]
                
                if product_index == 0 or pattern_idx == [] or last  == product_index - 1: 
                    pattern.append(atomic_num.GetSmarts())
                    pattern_idx.append(product_index)


                if atom == increment or atom == product_atoms_without_match[-1]:

                    idxs_product.append(pattern_idx)
                    patterns_product.append(pattern)

                
        atom_groups_product['product'] = Chem.MolToSmarts(product)
        atom_groups_product['patterns'] = patterns_product
        atom_groups_product['idxs'] = idxs_product


    for pattern in atom_groups_product['patterns']:
        
        pattern_with_match = []
        pattern_without_match = []
        




        # for reactant in reactants:
            
        #     substructure_atoms = [Chem.Atom(atom) for atom in pattern]
        #     substructure_smarts = Chem.MolToSmarts(Chem.MolFromAtoms(substructure_atoms))
        #     substructure = Chem.MolFromSmarts(substructure_smarts)
        #     reactant_substructure_pair = [reactant,substructure]

        #     if reactant.HasSubstructMatch(substructure):
        #         pattern_with_match['reactant']= reactant
        #         pattern_with_match['substructure']= substructure
        #     else:

             
        #     result =rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny)
        #     result_smarts = result.smartsString            

        # for increment in product_atom_increments:
        #     for atom in product_atoms_without_match:
                
        #         if product_atom_increments[product_atom_increments.index(increment) - 1] < atom <= increment:
        #             product_index = product_atom_increments.index(increment)


    results_dict['reactant atom count'] = reactant_atom_count
    results_dict['product atom count']= product_atom_count
    results_dict['match reactant atom count']= match_reactant_atom_count
    results_dict['atoms_without_match']= atoms_without_match
    results_dict['atoms_without_match_idx']= atoms_without_match_idx
    results_dict['reactant_atoms_with_match']= reactant_atoms_with_match
    results_dict['missing_patterns']= missing_patterns
    results_dict['product_atoms_with_match']= product_atoms_with_match
    results_dict['product_atoms_with_match_atom']= product_atoms_with_match_atom

    results_dict['product_atoms_without_match']= product_atoms_without_match
 
    results_dict['product_atom_increments']= product_atom_increments
    results_dict['atom_groups_product']= atom_groups_product
    
    print(results_dict) 
#IS SUBSTRUCTURE MATCH NOT COUNTING ATOM MATCH IF BOND IS different



# smiles_string = "CCO"
# reactant = 'CC(=O)O'
# Chem.MolFromSmarts(reactant)
# product_list = ['CC(=O)OC', 'CCO']
# for product in product_list:
#     Chem.MolFromSmarts(product)
agents = ['[H+]', '[Na+]']
rxn_smarts= "[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3]1=[C:4]"


reactant1 = Chem.MolFromSmiles('CC=O')
reactant2 = Chem.MolFromSmiles('CC(=O)CC')
product1 = Chem.MolFromSmiles('CCOC(=O)C')
product2 = Chem.MolFromSmiles('CCOCC(=O)CC')

reactants = [reactant1, reactant2]
reactant = reactant1
products = [product1, product2]
mols = [reactant1, product1]
find_closest_match(reactant, products)
find_best_combination(reactants, products)
# find_best_combination(reactants, products)



# Print the best combination
# for match_info in best_combination:
#     print(f"Reactant: {match_info['reactant']}, Product: {match_info['product']}")


# result= find_non_overlapping_matches(reactants, products)
# find_best_combination(matches_dict)
# print(result)
# Print the results
# for reactant, matches_info in result.items():
#     print(f"Reactant: {Chem.MolToSmiles(reactant)}")
#     for match_info in matches_info:
#         print(f"  Matched Product: {Chem.MolToSmiles(match_info['product'])}")
#         print(f"  Matching Substructure: {match_info['matching_substructure']}")
#         # print(f"  Reactant Atom Indices: {match_info['reactant_atom_indices']}")
#         # print(f"  Product Atom Indices: {match_info['product_atom_indices']}")
#         print()



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
#    [C:1](=[O:2])-[H:3]>>[C:1](-[H:3])-[O-:2]

# 'CCO.CC(=O)O>>CC(=O)OC.CCO'

rxn_dict = get_mols_from_rxn(rxn_smarts)
smarts= get_rxn_from_mols(rxn_dict)
# mol = get_indexed_smarts(smiles_string)
# highlight_similarity(rxn_dict)
# find_closest_match(reactant, product_list)
# print(result)
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