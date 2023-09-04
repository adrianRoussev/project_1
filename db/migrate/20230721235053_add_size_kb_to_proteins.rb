class AddSizeKbToProteins < ActiveRecord::Migration[7.0]
  def change
    add_column :proteins, :size_kb, :integer
  end
end
