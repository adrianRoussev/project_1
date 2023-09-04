class Intermediate < ApplicationRecord
    has_one :compoundable, as: :compoundable
    has_many :reaction_schemes, through: :compoundables, source: :compound, source_type: 'Compound'
  end