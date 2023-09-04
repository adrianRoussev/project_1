class CreateReactionSchemes < ActiveRecord::Migration[7.0]
  def change
    create_table :reaction_schemes do |t|
      t.integer :number_of_reactions
      t.integer :number_of_reatants
      t.string :similar_reaction_schemes
      t.timestamps
    end
  end
end
