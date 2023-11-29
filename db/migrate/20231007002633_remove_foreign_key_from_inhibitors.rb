class RemoveForeignKeyFromInhibitors < ActiveRecord::Migration[7.0]
  def change
    remove_foreign_key :inhibitors, :compounds
  end
end
