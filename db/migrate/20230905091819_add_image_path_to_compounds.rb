class AddImagePathToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :image_path, :string
  end
end
