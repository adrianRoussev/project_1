class Product < ApplicationRecord
  include Compoundable
    has_many :reactions
  end
  