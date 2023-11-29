class AddCompoundableToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_reference :compounds, :compoundable, polymorphic: true, index: true
  end
end
