class Reactant < ApplicationRecord
  include Compoundable
      has_many :reactions
  end