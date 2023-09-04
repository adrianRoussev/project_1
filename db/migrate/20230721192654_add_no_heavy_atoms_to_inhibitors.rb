class AddNoHeavyAtomsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :no_heavy_atoms, :integer
  end
end
