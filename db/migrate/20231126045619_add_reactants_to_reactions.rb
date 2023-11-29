class AddReactantsToReactions < ActiveRecord::Migration[7.0]
  def change
    add_column :reactions, :reactants, :string, array: true, default: []
    remove_column :reactions, :reactions, :string, array: true, default: []

  end
end
