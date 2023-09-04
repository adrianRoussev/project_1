class CreateCompounds < ActiveRecord::Migration[7.0]
  def change
    create_table :compounds do |t|
      t.float :molecular_weight
      t.string :molecular_formula
      t.string :smiles
      t.string :nucleophiles
      t.string :electrophiles
      t.integer :charge
      t.float :log_p 
      t.float :mass_available 
      t.string :name
      t.timestamps
    end
  end
end
