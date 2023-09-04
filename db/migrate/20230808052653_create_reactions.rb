class CreateReactions < ActiveRecord::Migration[7.0]
  def change
    create_table :reactions do |t|
      t.integer :number_of_products
      t.integer :number_of_reatants
      t.string :solvent
      t.string :catalysts
      t.integer :step_number
      t.string :type_of_reaction
      t.float :temperature_c
      t.string :similar_reaction_schemes

      t.timestamps
    end
  end
end
