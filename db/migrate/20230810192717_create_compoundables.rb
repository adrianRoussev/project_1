class CreateCompoundables < ActiveRecord::Migration[7.0]
  def change
    create_table :compoundables do |t|
      t.references :compound, polymorphic: true, null: false

      t.timestamps
    end
  end
end
