class AddCatalystArrayToReactions < ActiveRecord::Migration[7.0]
  def change
    remove_column :reactions, :catalysts, :string
    add_column :reactions, :catalysts, :string, array: true, default: []
  end
end
