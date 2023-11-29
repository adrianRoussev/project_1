class Compound < ApplicationRecord
belongs_to :compoundable, :polymorphic => true
belongs_to :reaction, optional: true
end

# #   compound = Compound.create(molecular_weight: 100, molecular_formula: 'C6H12O6')

# # # Creating an inhibitor associated with a compound
# # inhibitor = Inhibitor.create(name: 'Some Inhibitor')
# # inhibitor.compoundable = compound
# # inhibitor.save

# # # Accessing compounds associated with an inhibitor
# # inhibitor.compoundable

# # # Accessing inhibitors associated with a compound
# # compound.inhibitors


