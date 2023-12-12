class CreateCatalysts < ActiveRecord::Migration[7.0]
  def change
    create_table :catalysts do |t|
      t.references :reaction, null: false, foreign_key: true
      t.references :compound, null: false, foreign_key: true

      t.timestamps
    end
  end
end
