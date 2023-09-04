class CreateIntermediates < ActiveRecord::Migration[7.0]
  def change
    create_table :intermediates do |t|
      t.integer :product_number
      t.integer :step_number
      
      t.timestamps
    end
  end
end
