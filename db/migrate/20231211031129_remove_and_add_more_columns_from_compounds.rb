class RemoveAndAddMoreColumnsFromCompounds < ActiveRecord::Migration[7.0]
  def change
    remove_column :compounds, :fr_quatN, :integer
    remove_column :compounds, :fr_phenol_noOrthoHbond , :integer
    remove_column :compounds, :fr_nho, :integer
    add_column :compounds, :fr_quat_n, :integer
    add_column :compounds, :fr_phenol_no_ortho_h_bond , :integer
    add_column :compounds, :fr_nh0, :integer
  end
end
