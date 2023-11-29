class AddCompoundableTypeToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :compoundable_type, :string

  end
end
