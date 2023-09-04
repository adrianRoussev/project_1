class AddAllosteryToProteins < ActiveRecord::Migration[7.0]
  def change
    add_column :proteins, :allostery, :boolean
  end
end
