class CreateReactants < ActiveRecord::Migration[7.0]
  def change
    create_table :reactants do |t|
    t.string :reactive_functional_group
    t.float :mass_used
    t.float :concentration
    t.string :hazards
    t.float :melting_point_c 
    t.float :boiling_point_c 
    t.float :mass_available
      t.timestamps
    end
  end
end
