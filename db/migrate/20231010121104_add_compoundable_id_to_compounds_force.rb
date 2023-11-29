class AddCompoundableIdToCompoundsForce < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :compoundable_id, :bigint
  end
end
