class Reaction < ApplicationRecord
    has_many :reactants, class_name: 'Reactant', dependent: :destroy
    has_many :products, class_name: 'Product', dependent: :destroy
    has_many :catalysts, class_name: 'Catalyst', dependent: :destroy

    has_many :compounds, through: :reactants
    has_many :compounds, through: :products
    has_many :compounds, through: :catalysts

    accepts_nested_attributes_for :reactants
    accepts_nested_attributes_for :products
    accepts_nested_attributes_for :catalysts
    
    belongs_to :inhibitor

    

    # UNITS = {
    #     KILO: ['kilogram', 'KILOGRAM', 'Kilogram', 'KG', 'kg', 'Kg', 'kiloliter', 'KILOLITER', 'Kiloliter', 'KL', 'kl', 'Kl'],
    #     HECTO: ['hectogram', 'Hectogram', 'HECTOGRAM', 'hg', 'Hg', 'HG', 'hectoliter', 'Hectoliter', 'HECTOLITER', 'hl', 'Hl', 'HL'],
    #     DECA: ['decagram', 'Decagram', 'DECAGRAM', 'dag', 'Dag', 'DAG', 'decaliter', 'Decaliter', 'DECALITER', 'dal', 'Dal', 'DAL'],
    #     GRAM_LITRE: ['gram', 'Gram', 'GRAM', 'g', 'G','liter', 'Liter', 'LITER', 'L', 'l', 'LITRE', 'litre', ' '],
    #     CENTI: ['centigram', 'Centigram', 'CENTIGRAM', 'cg', 'Cg', 'CG', 'centiliter', 'Centiliter', 'CENTILITER', 'cl', 'Cl', 'CL'],
    #     MILLI: ['milligram', 'Milligram', 'MILLIGRAM', 'mg', 'Mg', 'MG', 'milliliter', 'Milliliter', 'MILLILITER', 'ml', 'Ml', 'ML'],
    #     MICRO: ['microgram', 'Microgram', 'MICROGRAM', 'ug', 'µg', 'Ug', 'UG', 'microliter', 'Microliter', 'MICROLITER', 'µl', 'Ul', 'UL'],
    #     NANO: ['nanogram', 'Nanogram', 'NANOGRAM', 'ng', 'Ng', 'NG', 'nanoliter', 'Nanoliter', 'NANOLITER', 'nl', 'Nl', 'NL'],
    #     PICO: ['picogram', 'Picogram', 'PICOGRAM', 'pg', 'Pg', 'PG', 'picoliter', 'Picoliter', 'PICOLITER', 'pl', 'Pl', 'PL'],
    #     FEMTO: ['femtogram', 'Femtogram', 'FEMTOGRAM', 'fg', 'Fg', 'FG', 'femtoliter', 'Femtoliter', 'FEMTOLITER', 'fl', 'Fl', 'FL']
    #     }


    # CONVERSION = {
    #     KILO: 10 ** 3,
    #     CENTA: 10 ** 2,
    #     DECA: 10 ** 1,
    #     GRAM_LITRE: 10 ** 0,
    #     DECI: 10 ** -1,
    #     CENTI: 10 ** -2,
    #     MILLI: 10 ** -3,
    #     MICRO: 10 ** -6,
    #     NANO: 10 ** -9,
    #     PICO: 10 ** -12,
    #     FEMTO: 10 ** -15
    #     }

    # SOLVENTS = {
    #     WATER: ['water', 'h2o', 'H2O', 'WATER', 'Water'],
    #     ACETONE: ['acetone', 'Acetone', 'acoac', 'AcO', 'ACO', 'aco', 'ACETONE', '(CH3)2CO', 'CH32CO', '(ch3)2co', 'ch32co'],
    #     METHANOL: ['methanol', 'Methanol', 'meoh', 'MeOh', 'MEOH', 'MeOH', 'METHANOL', 'ch3oh', 'CH3OH', '(CH3)OH',  '(ch3)oh'],
    #     ETHYL_ACETATE: ['ethyl acetate', 'Ethyl Acetate', 'EtOAc', 'ETOAC', 'etOAc', 'ETYLAC', 'ethylacetate', 'ETHYLACETATE'],
    #     DCM: ['dcm', 'DCM', 'dichloromethane', 'Dichloromethane', 'DICHLOROMETHANE', '(CH2)Cl2', 'CH2CL2', 'ch2cl2', '(ch2)cl2', 'DiChloroMethane', 'Methylene Dichloride', 'Methylene Chloride', 'methylene dichloride'],
    #     ETHANOL: ['ethanol', 'Ethanol', 'EtOH', 'ETOH', 'etoh', 'ch3ch2oh', 'CH3CH2OH', 'CH3-CH2-OH'],
    #     ETHER: ['ether', 'Ether', 'ETHER', 'diethyl ether', 'DIETHYL ETHER', 'Diethyl Ether', '(CH3CH2)2O', '(C2H5)2O', '(ch3ch2)2o', '(c2h5)2o'],
    #     CHLOROFORM: ['chloroform', 'Chloroform', 'CHCl3', 'chcl3', 'CHCL3'],
    #     TOLUENE: ['toluene', 'Toluene', 'TOLUENE'],
    #     HEXANE: ['hexane', 'Hexane', 'HEXANE']
    #         }    
  
    # DENSITIES = {
    #     SOLVENTS[:WATER].each { |water| water => 1.00 },
    #     SOLVENTS[:ACETONE].each { |acetone| acetone => 0.79 },
    #     SOLVENTS[:METHANOL].each { |methanol| methanol => 0.79 },
    #     SOLVENTS[:ETHYL_ACETATE].each { |ethyl_acetate| ethyl_acetate => 0.902 },
    #     SOLVENTS[:DCM].each { |dcm| dcm => 1.33 },
    #     SOLVENTS[:ETHANOL].each { |ethanol| ethanol => 0.789 },
    #     SOLVENTS[:ETHER].each { |ether| ether => 0.713 },
    #     SOLVENTS[:CHLOROFORM].each { |chloroform| chloroform => 1.49 },
    #     SOLVENTS[:TOLUENE].each { |toluene|toluene => 0.867 },
    #     SOLVENTS[:HEXANE].each { |hexane| hexane => 0.659 }
    #         } 
    
    
        
    def create_reaction(products, reactants, catalysts)

        products_array = get_smiles_array(products)
        catalysts_array = get_smiles_array(catalysts)
        reactants_array = get_smiles_array(reactants)

    
        products_str = join_smiles(products_array)
        catalysts_str = join_smiles(catalysts_array)
        reactants_str = join_smiles(reactants_array)

        reaction = "#{reactants_str}>#{catalysts_str}>#{products_str}"

        reaction
    end


    def get_smiles_array(compound_array)
        smiles_array = []

        compound_array.each do |x|
            smiles = x.compound.smiles
            smiles_array << smiles
        end

        smiles_array
    end

    def create_reaction_dict(products, reactants, catalysts)
        products_array = get_smiles_array(products)
        catalysts_array = get_smiles_array(catalysts)
        reactants_array = get_smiles_array(reactants)
        rxn_dict= {
                        "reactants"=>  reactants_array,
                        "products"=> products_array,
                        "catalysts"=> catalysts_array
                    }
       
        rxn_dict

    end

    def get_compound_quantity(compound_array)
        compound_quantity = {}

        compound_array.each do |x| 
            if x.is_a?(Reactant) || x.is_a?(Catalyst)
                quantity = x.mass_used
            elsif x.is_a?(Product)
                quantity = x.actual_yeild
            end   
            compound_quantity[x.compound.smiles]  = quantity
        end

        compound_quantity
    end

    def convert_mass_or_vol(value, initial_unit, final_unit)
    
        conversion = {'kilo': 10 ** (3),
                          'centa': 10 ** (2),
                          'deca': 10 ** (1),
                          'gram': 10 ** (0),
                          'deci': 10 ** (-1),
                          'centi': 10 ** (-2),
                          'milli': 10 ** (-3),
                          'micro': 10 ** (-6),
                          'nano': 10 ** (-9),
                          'pico': 10 ** (-12),
                        }
        converted_value = (value/conversion[initial_unit])* conversion[final_unit]
        converted_value
    
    end 
    
    def calculate_n(mass, multiplier, power_number) #n = numerator... should probably rename so doesnt get mixed up with mols

        value = mass * (multiplier ** power_number)

        value
        
    end

    def calculate_mols(molecular_weight, mass_in_g) #should this be in compound not reaction?

        mols = mass_in_g/molecular_weight

        mols
    end


    def mass_from_mols(mols, molecular_weight)

        mass = mols * molecular_weight

        mass
    end 

    def solvent_weight(reaction, weight, initial_weight_units)



        density_units = "g/ml"
        density = densities[reaction.solvent.to_s.downcase]
        
        correct_weight_units = "gram"
        corrected_weight = convert_mass_or_vol(weight, intitial_weight_units, correct_weight_units)
        volume = weight/density

        volume

    end
        

    def concentration_formula(compound, reaction, power_number, type)

        numerators = {'mass/vol': compound.mass,
            'Molarity':  calculate_mols(compound.molecular_weight, compound.mass),
            'mass %': calculate_n(compound.mass, 100, 1) , 
            'volume %': calculate_n(reaction.solvent_volume, 100, 1) ,
            'ppt/m/b': calculate_n(compound.mass, 10, power_number),
            'normality': calculate_n(compound.mass, solvent_weight(reaction, weight), -1)  
            }
    
        denominators = {'mass/vol': vol,
            'Molarity': reaction.solvent_volume,
            'mass %': (solvent_weight(reaction, weight)+ compound.mass),
            'volume %': reaction.solvent_volume ,
            'ppt/m/b': (solvent_weight(reaction, weight)+ compound.mass),
            'normality': reaction.solvent_volume
            }

        n = numerators[type]
        d = denominators[type]

        numerator = convert_mass_or_vol(n, initial_unit_n, final_unit_n)
        denominator = convert_mass_or_vol(d, initial_unit_d, final_unit_d)
        

    
        concentration_value = numerator/denominator
        units = "#{final_unit_n}/#{final_unit_d}"

        concentration = "#{concentration_value} #{units}"
    
        concentration

    end

    

    def find_least_mols( reactants) #finds single eq when mols or mass is known for every reactant 

        reactant_mols = get_reactant_mols_hash(reactants)
    
        sorted_reactant_mols = reactant_mols.sort_by{|key, value| value}
        least_mols = sorted_reactant_mols.values.first
        
        least_mols
        end

    def eq_for_compound(reactants, compound_of_interest)

        least_mols = find_least_mols(reactants)
    
        equivalents = fdiv(compound_of_interest.mols,least_mols)

        equivalents
    end

    def  mol_eq_from_mol_and_eq(reactants) #use when eq and mass/mols are known

        single_eq_mols_array = []
        reactants.each do |reactant|

            if reactant.mols && reactant.equivalents &&reactants_with_mass.empty? || reactant.mass && reactant.equivalents &&reactants_with_mass.empty?
                single_eq_mols =  fdiv(mols,reactant.equivalents)
                single_eq_mols_array << single_eq_mols
           
            end
        end
        
        single_eq_mol = reactants_with_mass[0]
        single_eq_mol
    end


    def find_single_eq_mol_from_unknown(reactants) # --> use when determining conditions are unknown 

        if check_reactants_mols_or_mass(reactants) == true
            single_eq_mols = find_least_mols(reactants)
        
        elsif check_for_mols_or_mass_and_eq(reactants) == true
            single_eq_mols = mol_eq_from_mol_and_eq(reactants)

        end
        
        single_eq_mols  
    end
        
        #i have an array of reactants for which each reactant has an attribute reactant.mol, reactant.mass and reactant.eq . these attributes are not always populated. i want to write a loop that checks wheather any reactants exist that have either both reactant.mols & reactant.eq or reactant.mass & reactant.eq being assigned. if yes stop loop, if no continue through list of reactants. if none of the reactants returned true then return false. if all reactants meet requirements ie. the loop is complete without interuption, return true 


    def mass_from_eq(reaction, compound_of_interest)

        single_eq_mol = mol_eq_from_mol_and_eq(reactants)
        mols = compound_of_interest.equivalents * single_eq_mol

        mass = mass_from_mols(mols, molecular_weight)

        mass

    end

    def find_max_yeild(product, reactants)
        
        limiting_mols = find_single_eq_mol_from_unknown(reactants) 

        #if product.composing_reactants.include?(limiting_reactant) && if !reactants.include?(product)
        mass = mass_from_mols(limiting_mols, molecular_weight)

        mass

        # ** need to consider stoichiometry ie. if reactant takes 2x limiting mol then limiting must be divided 
    
    end 



        

        





        




        
        
    # def get_compound_data

    # def get_mass_from_mol

    # def  get_mass_from_eq

    # def get_mass_from_concentration 

    # def get_max_mols_from_eq

    
   

    private 

    def join_smiles(smiles_array)
        compound_str = smiles_array.join(".")
        compound_str
    end

    def vol_from_density(mass, density)
        volume = mass/density
        volume
    end



  
  def convert_mass_or_vol(value, initial_unit, final_unit)
    converted_value = (value / CONVERSION_RATES[initial_unit]) * CONVERSION_RATES[final_unit]
    converted_value
  end
  
  # Solvent weight calculation method
  def solvent_weight(reaction, weight, initial_weight_units)
    density = DENSITIES[reaction.solvent.downcase.to_sym]
    correct_weight_units = Unit::G
    corrected_weight = convert_mass_or_vol(weight, initial_weight_units, correct_weight_units)
    volume = weight / density
    volume
  end
  
  def get_reactant_mols_hash(reactants) # finds mols of reactants from if mass or mols is present

    reactant_mols = {}
    
    reactants.each do |reactant| 
        if reactant.mols?
            mols = reactant.mols
        else
            mols = reactant.calculate_mols(reactant.compound.molecular_weight, reactant.mass)
        end
        
        reactant_mols[reactant] = mols
    end

    reactant_mols
  end

  def check_for_mols_or_mass(reactants)
    reactants.each do |reactant|
      unless reactant[:mol] || reactant[:mass]
        return false 
      end
    end
    true 
  end

  def check_for_mols_or_mass_and_eq(reactants)

    reactants.each do |reactant|
      if (reactant[:mol] && reactant[:eq]) || (reactant[:mass] && reactant[:eq])
        return true 
      end
    end
    false 
  end
end