from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS
from itertools import permutations
from itertools import product
from itertools import combinations



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


def find_value_combination(input_set, input_list):
    # print("input", input_set, input_list)
    total_matches = []
    for item in input_list:
        for pattern, matches in item.items():
            for match in matches:
                total_matches.append(match)
    # print("total_matches:", total_matches)
    
    count = len(total_matches)

    index_set = set(range(1,count))

    
    match_lengths = []
    for each in total_matches:
        match_lengths.append(len(each))

    # print(match_lengths)

    min_length = min(match_lengths)
    max_length = max(match_lengths)
    product_length = len(input_set)

    length_by_min = product_length//min_length
    length_by_max = product_length//max_length

    combination_n_elements = list(range(length_by_max, length_by_min + 1))

    # print(combination_n_elements)

    combinations_dict = {}
    
    for n_element in combination_n_elements:
        
        unique_pairs = combinations(index_set, n_element)
        combinations_dict[n_element] = unique_pairs

    # print(unique_pairs)
    unique_sets = []
    for each in combinations_dict.values():
        for group in each:
            unique_sets.append(group)
    
    # print(unique_sets)
    
    list_combinations = []
    for pair in unique_sets:
        list_pair = []
        for index in pair:
            list_pair.append(total_matches[index])
        # if list_combinations:   
        #     for element in list_combinations:
        #         set_element = {tuple(sorted(sub)) for sub in element}
        #         if all(set_element != {tuple(sorted(sub)) for sub in sublist} for sublist in list_combinations):
        #             list_combinations.append(list_pair)
        # else:
        list_combinations.append(list_pair)




    # print('list combination:', list_combinations)
    
    value_list = []
    for combination in list_combinations:
        value_set = []
        for number_list in combination:
            for value in number_list:
                value_set.append(value)
        value_list.append(value_set)
    
    missing_values_list = []
    for value_set in value_list:
        missing_values = []
        for value in input_set:
            if value not in value_set:
                missing_values.append(value)
        missing_values_list.append(missing_values)

    # print('missing values list:', missing_values_list)
    
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
    # print(missing_values_dict)
   
    # for each in missing_values_dict.keys():
    #     keys_list.append(each)
    # if keys_list: #remove after 'if no_missing_values_index == []:' indenting is fixed
    #     min_key = min(keys_list)

    #     min_dict = {key: value for key, value in missing_values_dict.items() if key == min_key}

        # print(min_dict)

        
        
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
    consecutive_pairs = []
    non_consecutive_pairs = []
    combinations_with_missing = {}
    if no_missing_values_index:
        for index in no_missing_values_index:
            selected_combinations.append(list_combinations[index])
    if missing_values_dict: #change to elif
       
        if len(missing_values_dict) ==1:
            selected_combinations.append(missing_values_dict.values())
        else:
            for each in missing_values_dict.values(): 
                unknown_list = each[0] 
                is_consecutive = True
                unknown_list.sort()  
                
                for i in range(len(unknown_list) - 1):
                    if unknown_list[i] + 1 != unknown_list[i + 1]:
                        is_consecutive = False
                        non_consecutive_pairs.append(unknown_list)
                        
                    break
                if is_consecutive == True:
                    consecutive_pairs.append(unknown_list)
                    
    combinations_with_missing['consecutive pairs'] = consecutive_pairs
    combinations_with_missing['non-consecutive pairs'] = non_consecutive_pairs

    results_dict = {}
    results_dict['selected_combinations'] = selected_combinations
    results_dict['combinations with missing values'] = combinations_with_missing

    # print(results_dict)
    return results_dict

                        
# - Prioritize Missing Atoms Matching Solvent or Catalyst: If any missing atoms match the solvent or catalyst, 
#   prioritize those reactant pairs. This ensures that environmental factors are considered in the reactant selection process.
#   (consecutive)

# - Prioritize Reactions with OH or C=O with Nucleophile Next to It: Reactant pairs that can potentially form OH or C=O groups with a 
#   nucleophile nearby should be prioritized. This accounts for common reaction environments where these functional groups are likely to occur. 
#   (single missing or 2 consecutive) 

# - Prioritize Single Missing Atom Where Single Atom is Not Carbon: If a single atom is missing and it's not carbon, prioritize reactant
#   pairs that can fill this gap. This helps to ensure completeness in the product formation. (single)

#- Prioritize Reactions with No Missing Carbons Unless O=C=O Could Have Formed: If no carbons are missing in the product, prioritize 
#   reactant pairs accordingly, unless the possibility of forming O=C=O exists. This accounts for potential variations in the reaction mechanisms.
#   (single or non-consecutive)
                        
    
# Button for 'if these predictions are wrong' search database for similar reactions and select most similar reactants

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


def find_best_combination(rxn_smarts):

    products = rxn_smarts.GetProducts()
    reactants = rxn_smarts.GetReactants()
    catalysts = rxn_smarts.GetAgents()
    

    matches_dict = {}
    product_combos= {}
    product_result = {}
    for reactant in reactants:
        
        matches_dict[Chem.MolToSmarts(reactant)] = find_closest_match(reactant, products)

    print('matches_dict', matches_dict)
    for product in products:
        product_combos_list = []
        count = 0
        product_atom_idxs = []
       
        for atom in product.GetAtoms():
            product_atom_idxs.append(count)
            count+=1
        # product_atom_idxs = list(range(product.GetAtoms))
       
        for reactant in reactants:
            combos = matches_dict[Chem.MolToSmarts(reactant)][Chem.MolToSmarts(product)]
            product_combos_list.append({Chem.MolToSmarts(reactant): combos})
        
        product_combos[Chem.MolToSmarts(product)] = product_combos_list
        # print("input params", product_atom_idxs, product_combos_list)
        print(product_combos)
        result = find_value_combination(product_atom_idxs, product_combos_list)
        product_result[Chem.MolToSmarts(product)] = result

    print('product_result:', product_result)
    print('product_combos:', product_combos)
    


    product_reactants = {}
    product_options = {}
    product_options_w_atomic_n = {}
    for product in products:
        if product_result[Chem.MolToSmarts(product)]['selected_combinations']:
           
            product_reactants[Chem.MolToSmarts(product)] = product_result[Chem.MolToSmarts(product)]['selected_combinations']


        # else: 
            missing_values_dict = product_result[Chem.MolToSmarts(product)]['combinations with missing values'] 
            reactants_for_product= []
            for each in product_combos[Chem.MolToSmarts(product)]:
                reactants_for_product.append(each.keys())          
    
            product_options[Chem.MolToSmarts(product)] = {'reactants': [reactants_for_product],'missing product atoms': missing_values_dict}
    
        consecutive_pairs = product_options[Chem.MolToSmarts(product)]['missing product atoms']['consecutive pairs']
        non_consecutive_pairs = product_options[Chem.MolToSmarts(product)]['missing product atoms']['non-consecutive pairs']
        reactants_for_product = product_options[Chem.MolToSmarts(product)]['reactants']

        consecutive_atom_groups = []
        non_consecutive_atom_groups = []

        for reactant_possibility in reactants_for_product:
            if consecutive_pairs:
    
                for pair in consecutive_pairs:
                    missing_indecies =[]
                    
                    for element in pair:
                        missing_indecies.append(product.GetAtomWithIdx(element).GetAtomicNum())
                    consecutive_atom_groups.append(missing_indecies)
            if non_consecutive_pairs:
                for group in non_consecutive_pairs:
                    missing_indecies =[]
                    
                    for element in group:
                        missing_indecies.append(product.GetAtomWithIdx(element).GetAtomicNum())
                    non_consecutive_atom_groups.append(missing_indecies)

        product_options_w_atomic_n[Chem.MolToSmarts(product)] = {'reactants':reactants_for_product, 'non_consecutive': non_consecutive_atom_groups, 'consecutive': consecutive_atom_groups}

        ordered_product_options = []
        
        for each in product_options_w_atomic_n[Chem.MolToSmarts(product)]['consecutive']:
            all_atoms_list = []
            matches = sum(1 for x in each if x in all_atoms_list)
            
            for catalyst in catalysts:

                catalyst_atoms = catalyst.GetAtoms()
                all_atoms_list = []
            
                for each in catalyst_atoms:
                    all_atoms_list.append(each.GetAtomicNum())

                    catalyst_percentage_match = (matches / ((len(each)+len(all_atoms_list))/2)) * 100

                if catalyst_percentage_match >= 80:
                    ordered_product_options.append(each)

            common_additions = [[8], [6,8,8], [9], [17], [35], [53]]
            
            for common_mol in common_additions:
                all_atoms_list = common_mol
                common_percentage_match = (matches / ((len(each)+len(all_atoms_list))/2)) * 100

                if common_percentage_match >= 80:
                    ordered_product_options.append(each)

        for each in product_options_w_atomic_n[Chem.MolToSmarts(product)]['non_consecutive']:

            common_additions = [[8], [9], [17], [35], [53]]

            for common_atom in common_additions:

                if common_atom in each and 6 not in each: 
                    ordered_product_options.append(each)

    print(product_reactants)
    print(ordered_product_options)
    
    if product_reactants:
        return product_reactants
    
    else: 
        return ordered_product_options
        



    



    # mol = 'C(=O)=O'
    
    # all_atoms = (Chem.MolFromSmarts(mol)).GetAtoms()
    # all_atoms_list = []
    # for each in all_atoms:
    #     all_atoms_list.append(each.GetAtomicNum())
    # print(all_atoms_list)



                

                


#   Prioritize Missing Atoms Matching Solvent or Catalyst: If any missing atoms match the solvent or catalyst, 
#   prioritize those reactant pairs. This ensures that environmental factors are considered in the reactant selection process.
#   (consecutive)

# - Prioritize Reactions with OH or C=O with Nucleophile Next to It: Reactant pairs that can potentially form OH or C=O groups with a 
#   nucleophile nearby should be prioritized. This accounts for common reaction environments where these functional groups are likely to occur. 
#   (single missing or 2 consecutive) 

# - Prioritize Single Missing Atom Where Single Atom is Not Carbon: If a single atom is missing and it's not carbon, prioritize reactant
#   pairs that can fill this gap. This helps to ensure completeness in the product formation. (single)

#- Prioritize Reactions with No Missing Carbons Unless O=C=O Could Have Formed: If no carbons are missing in the product, prioritize 
#   reactant pairs accordingly, unless the possibility of forming O=C=O exists. This accounts for potential variations in the reaction mechanisms.
#   (single or non-consecutive)
                        
    
# Button for 'if these predictions are wrong' search database for similar reactions and select most similar reactants












#EXAMPLES


# smiles_string = "CCO"
# reactant = 'CC(=O)O'
# Chem.MolFromSmarts(reactant)
# product_list = ['CC(=O)OC', 'CCO']
# for product in product_list:
#     Chem.MolFromSmarts(product)
# agents = ['[H+]', '[Na+]']
# rxn_smarts_1= "[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3]1=[C:4]"
rxn_smarts="[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3]1=[C:4]"
rxn =Chem.rdChemReactions.ReactionFromSmarts(rxn_smarts)



reactant1 = 'CC=O'
reactant2 = 'CC(=O)CC'
product1 = 'CCOC(=O)C'
product2 = 'CCOCC(=O)CC'

reactants = [reactant1, reactant2]
reactant = reactant1
products = [product1, product2]
mols = [reactant1, product1]
catalysts = []
rxn_dict = {
        'reactants': reactants,
        'products': products,
        'catalysts': catalysts
    }
find_best_combination(rxn)



# # Print the best combination
# # for match_info in best_combination:
# #     print(f"Reactant: {match_info['reactant']}, Product: {match_info['product']}")


# #

# # 1. Aldol Condensation:
# #    [C:1](=[O:2])-[CH2:3]-[C:4](=[O:5])>>[C:1](=[O:2])-[C:3]([OH:6])-[C:4](=[O:5])

# # 2. Grignard Reaction:
# #    [Mg:1].[C:2](=[O:3])>>[C:2]([Mg:1])[OH2:3]
# #[Mg][CX3:1][X:2].[C:3](=[O:4])[H:5]>>[C:3][X:2][CX3:1][O:4][H:5]
# # 3. Ring Closure - Diels-Alder:
# #    [C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3]1=[C:4]

# # 4. Substitution - SN2:
# #    [C:1]-[Br:2].[N:3]>>[C:1]-[N:3].[Br:2]

# # 5. Oxidation of Alcohol:
# #    [C:1]-[OH:2]>>[C:1](=[O:3])-[OH:2]

# # 6. Reduction of Ketone:
# #    [C:1](=[O:2])-[H:3]>>[C:1](-[H:3])-[O-:2]

# # 'CCO.CC(=O)O>>CC(=O)OC.CCO'

# rxn_dict = get_mols_from_rxn(rxn_smarts)
# smarts= get_rxn_from_mols(rxn_dict)
# print("Reaction SMARTS:", smarts)

