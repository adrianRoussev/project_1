class AddPolymorphicAssociationToCompoundables < ActiveRecord::Migration[7.0]
  def change
    add_column :compoundables, :compoundable_type, :string
    add_column :compoundables, :compoundable_id, :bigint
    add_index :compoundables, [:compoundable_type, :compoundable_id], name: 'index_compoundables_on_compoundable'
  end
end
