class Reaction < ApplicationRecord
    has_many :reactants, class_name: 'Reactant'
    has_many :products, class_name: 'Product'
    has_many :catalysts, class_name: 'Catalyst'

    has_many :compounds, through: :reactants
    has_many :compounds, through: :products
    has_many :compounds, through: :catalysts

    accepts_nested_attributes_for :reactants
    accepts_nested_attributes_for :products
    accepts_nested_attributes_for :catalysts
    
    belongs_to :inhibitor

    def create_reaction(products, reactants, catalysts)

        products_array = get_smiles_array(products)
        catalysts_array = get_smiles_array(catalysts)
        reactants_array = get_smiles_array(reactants)

       
        products_str = join_smiles(products_array)
        catalysts_str = join_smiles(catalysts_array)
        reactants_str = join_smiles(reactants_array)

        reaction = "#{products_str}>#{catalysts_str}>#{reactants_str}"

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


end

