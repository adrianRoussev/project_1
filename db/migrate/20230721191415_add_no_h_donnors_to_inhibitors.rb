class AddNoHDonnorsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :no_h_donnors, :integer
  end
end
