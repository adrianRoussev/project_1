class AddReactionsAndCompoundsToReactants < ActiveRecord::Migration[7.0]
  def change
    add_reference :reactants, :reaction, null: false, foreign_key: true
    add_reference :reactants, :compound, null: false, foreign_key: true

  end
end
