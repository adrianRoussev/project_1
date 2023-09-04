class AddOncogenicToProteins < ActiveRecord::Migration[7.0]
  def change
    add_column :proteins, :oncogenic, :boolean
  end
end
