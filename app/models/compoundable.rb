class Compoundable < ApplicationRecord
    belongs_to :compound
    belongs_to :compoundable, polymorphic: true
  end