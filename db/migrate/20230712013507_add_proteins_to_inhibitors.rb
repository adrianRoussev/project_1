class AddProteinsToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_reference :inhibitors, :protein, null: false, foreign_key: true
  end
end
