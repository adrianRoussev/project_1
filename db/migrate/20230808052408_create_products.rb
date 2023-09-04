class CreateProducts < ActiveRecord::Migration[7.0]
  def change
    create_table :products do |t|
    t.string :hazards
    t.float :melting_point_c
    t.float :boiling_point_c
    t.string :reported_reactions
    t.integer :reported_yeild
    t.integer :actual_yeild
    t.float :obtained_mass

      t.timestamps
    end
  end
end
