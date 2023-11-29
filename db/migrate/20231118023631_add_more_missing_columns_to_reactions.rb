class AddMoreMissingColumnsToReactions < ActiveRecord::Migration[7.0]
  def change
    add_column :reactions, :reactions, :string, array: true, default: []
    add_column :reactions, :products, :string, array: true, default: []
    add_column :reactions, :generated_products, :string, array: true, default: []
  end
end
