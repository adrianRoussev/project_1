class AddCompoundsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_reference :inhibitors, :compound, null: false, foreign_key: true
  end
end
