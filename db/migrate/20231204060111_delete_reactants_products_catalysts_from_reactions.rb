class DeleteReactantsProductsCatalystsFromReactions < ActiveRecord::Migration[7.0]
  def change
    remove_column :reactions, :reactants, :string, array: true, default: []
    remove_column :reactions, :products, :string, array: true, default: []
    remove_column :reactions, :catalysts, :string, array: true, default: []
  end
end
