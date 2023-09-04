# db/migrate/xxxxxx_add_compound_properties_to_inhibitors.rb
class AddNewColumnToInhibitors < ActiveRecord::Migration[7.0]
  def change
    add_column :inhibitors, :molecular_formula, :string
    add_column :inhibitors, :molecular_weight, :string
    add_column :inhibitors, :canonical_smiles, :string
    add_column :inhibitors, :isomeric_smiles, :string
    add_column :inhibitors, :inchi, :string
    add_column :inhibitors, :inchi_key, :string
    add_column :inhibitors, :iupac_name, :string
    add_column :inhibitors, :title, :string
    add_column :inhibitors, :x_log_p, :string
    add_column :inhibitors, :exact_mass, :string
    add_column :inhibitors, :monoisotopic_mass, :string
    add_column :inhibitors, :tpsa, :string
    add_column :inhibitors, :complexity, :string
    add_column :inhibitors, :charge, :string
    add_column :inhibitors, :h_bond_donor_count, :string
    add_column :inhibitors, :h_bond_acceptor_count, :string
    add_column :inhibitors, :rotatable_bond_count, :string
    add_column :inhibitors, :heavy_atom_count, :string
    add_column :inhibitors, :isotope_atom_count, :string
    add_column :inhibitors, :atom_stereo_count, :string
    add_column :inhibitors, :defined_atom_stereo_count, :string
    add_column :inhibitors, :undefined_atom_stereo_count, :string
    add_column :inhibitors, :bond_stereo_count, :string
    add_column :inhibitors, :defined_bond_stereo_count, :string
    add_column :inhibitors, :undefined_bond_stereo_count, :string
    add_column :inhibitors, :covalent_unit_count, :string
    add_column :inhibitors, :patent_count, :string
    add_column :inhibitors, :patent_family_count, :string
    add_column :inhibitors, :literature_count, :string
  end
end
