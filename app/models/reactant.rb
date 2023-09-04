class Reactant < ApplicationRecord
    has_one :compoundable, as: :compoundable
    has_many :reactions
  end