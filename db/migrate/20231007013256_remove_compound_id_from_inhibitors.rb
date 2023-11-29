class RemoveCompoundIdFromInhibitors < ActiveRecord::Migration[7.0]
  def change
    remove_column :inhibitors, :compound_id, :bigint, null: false
  end
end

