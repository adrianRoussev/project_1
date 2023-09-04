class Reaction < ApplicationRecord
    has_many :products, through: :compounds
    has_many :reactants, through: :compounds

    def join_compounds(*smiles)
      rxn_compounds = smiles.join(".")
      rxn_compounds
    end

    def create_reaction(products, reactants, catalysts)
      reaction = "#{products_str}>#{catalysts}>#{reactants_str}"
      reaction
    end
end


