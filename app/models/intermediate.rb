class Intermediate < ApplicationRecord
  include Compoundable
    # has_many :reaction_schemes, through: :compoundables, source: :compound, source_type: 'Compound'
  end