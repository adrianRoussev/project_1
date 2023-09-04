class AddNoImagePathToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :image_path, :string
  end
end
