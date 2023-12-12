class AddMolecularParamsToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :molecular_params, :string
  end
end
