class AddMissingColumnsToReactions < ActiveRecord::Migration[7.0]
  def change
    add_column :reactions, :target_product_name, :string
    add_column :reactions, :reaction_smiles, :string


  end
end
