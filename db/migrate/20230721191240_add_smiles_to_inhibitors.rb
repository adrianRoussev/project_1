class AddSmilesToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :smiles, :string
  end
end
