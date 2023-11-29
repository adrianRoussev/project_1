class RemoveColumnsFromCompounds < ActiveRecord::Migration[7.0]
  def change
    remove_column :compounds, :compoundable_type, :string
    remove_column :compounds, :compoundable_id, :string

  end
end
