class ReactionScheme < ApplicationRecord
    belongs_to :inhibitor
    has_many :reactions
    has_many :intermediates, through: :reactions
  end

  def create_reaction_scheme (*reaction)
      reaction_scheme_array = []
      reaction_scheme_array << reaction
      reaction_scheme
  end 

  