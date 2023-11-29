class AddInhibitorsToReactions < ActiveRecord::Migration[7.0]
  def change
    add_reference :reactions, :inhibitor, null: false, foreign_key: true

  end
end

 