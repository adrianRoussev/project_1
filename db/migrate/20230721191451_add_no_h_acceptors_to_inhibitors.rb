class AddNoHAcceptorsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :no_h_acceptors, :integer
  end
end
