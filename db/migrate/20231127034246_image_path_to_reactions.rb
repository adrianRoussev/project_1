class ImagePathToReactions < ActiveRecord::Migration[7.0]
  def change
    add_column :reactions, :image_path, :string
  end
end
