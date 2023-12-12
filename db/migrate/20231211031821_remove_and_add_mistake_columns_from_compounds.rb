class RemoveAndAddMistakeColumnsFromCompounds < ActiveRecord::Migration[7.0]
  def change
    remove_column :compounds, :fr_quat_n, :integer
    remove_column :compounds, :fr_phenol_no_ortho_h_bond , :integer
    add_column :compounds, :fr_quatn, :integer
    add_column :compounds, :fr_phenol_noorthohbond , :integer
  end
end
