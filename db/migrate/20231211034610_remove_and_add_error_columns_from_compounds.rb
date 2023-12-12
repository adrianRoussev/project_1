class RemoveAndAddErrorColumnsFromCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :bcut2d_mrlow , :float
    remove_column :compounds, :bcut2d_mrlo , :float
  end
end
