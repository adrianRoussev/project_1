class AddCidToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :cid, :string
  end
end
