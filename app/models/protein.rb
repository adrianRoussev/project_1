class Protein < ApplicationRecord
    has_many :compoundables, as: :compoundable
  end
  