class RemoveSimilarReactionSchemesFromReactions < ActiveRecord::Migration[7.0]
  def change
    remove_column :reactions, :similar_reaction_schemes, :string
  end
end
