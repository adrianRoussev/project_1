class Reaction < ApplicationRecord
    has_many :compounds
    belongs_to :inhibitor

    def create_reaction(products, reactants, catalysts)
        products_str = products.
      reaction = "#{products_str}>#{catalysts_str}>#{reactants_str}"
      reaction
    end

private 

def join_smiles(*smiles)
   compound_str = smiles.join(".")
   compound_str
  end

end

