class AddTargetSiteToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :target_site, :string
  end
end
