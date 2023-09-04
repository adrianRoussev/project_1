class Compound < ApplicationRecord
    has_many :compoundables
    has_many :inhibitors, through: :compoundables, source: :compoundable, source_type: 'Inhibitor'
    has_many :intermediates, through: :compoundables, source: :compoundable, source_type: 'Intermediate'
    has_many :products, through: :compoundables, source: :compoundable, source_type: 'Product'
    has_many :reactants, through: :compoundables, source: :compoundable, source_type: 'Reactant'
    
    has_many :reactions, through: :products
    has_many :reactions, through: :reactants
 


#   compound = Compound.create(molecular_weight: 100, molecular_formula: 'C6H12O6')

# # Creating an inhibitor associated with a compound
# inhibitor = Inhibitor.create(name: 'Some Inhibitor')
# inhibitor.compoundable = compound
# inhibitor.save

# # Accessing compounds associated with an inhibitor
# inhibitor.compoundable

# # Accessing inhibitors associated with a compound
# compound.inhibitors


end