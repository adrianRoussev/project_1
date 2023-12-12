class AddReactionsAndCompoundsToProducts < ActiveRecord::Migration[7.0]
  def change
    add_reference :products, :reaction, null: false, foreign_key: true
    add_reference :products, :compound, null: false, foreign_key: true

  end
end
