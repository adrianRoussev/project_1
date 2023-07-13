class CreateInhibitor < ActiveRecord::Migration[7.0]
  def change
    create_table :inhibitors do |t|
      t.string :name
      t.string :structure
      t.boolean :lipinskis

      t.timestamps
    end
  end
end
