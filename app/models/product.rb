class Product < ApplicationRecord
  include Compoundable
    belongs_to :reaction
    belongs_to :compound
end
  