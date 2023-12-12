class AddMassUsedToCatalysts < ActiveRecord::Migration[7.0]
  def change
    add_column :catalysts, :mass_used, :float
  end
end
