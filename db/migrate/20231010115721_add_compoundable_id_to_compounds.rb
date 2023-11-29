class AddCompoundableIdToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :compoundable_id, :string

  end
end
