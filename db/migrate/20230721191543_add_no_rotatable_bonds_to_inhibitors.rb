class AddNoRotatableBondsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :no_rotatable_bonds, :integer
  end
end
