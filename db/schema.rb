# This file is auto-generated from the current state of the database. Instead
# of editing this file, please use the migrations feature of Active Record to
# incrementally modify your database, and then regenerate this schema definition.
#
# This file is the source Rails uses to define your schema when running `bin/rails
# db:schema:load`. When creating a new database, `bin/rails db:schema:load` tends to
# be faster and is potentially less error prone than running all of your
# migrations from scratch. Old migrations may fail to apply correctly if those
# migrations use external dependencies or application code.
#
# It's strongly recommended that you check this file into your version control system.

ActiveRecord::Schema[7.0].define(version: 2023_08_10_192717) do
  # These are extensions that must be enabled in order to support this database
  enable_extension "plpgsql"

  create_table "compoundables", force: :cascade do |t|
    t.string "compound_type", null: false
    t.bigint "compound_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.index ["compound_type", "compound_id"], name: "index_compoundables_on_compound"
  end

  create_table "compounds", force: :cascade do |t|
    t.float "molecular_weight"
    t.string "molecular_formula"
    t.string "smiles"
    t.string "nucleophiles"
    t.string "electrophiles"
    t.integer "charge"
    t.float "log_p"
    t.float "mass_available"
    t.string "name"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "inhibitors", force: :cascade do |t|
    t.string "name"
    t.string "structure"
    t.boolean "lipinskis"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.bigint "protein_id", null: false
    t.string "target_site"
    t.string "smiles"
    t.integer "no_h_donnors"
    t.integer "no_h_acceptors"
    t.integer "no_rotatable_bonds"
    t.integer "no_heavy_atoms"
    t.string "image_path"
    t.string "molecular_formula"
    t.string "molecular_weight"
    t.string "canonical_smiles"
    t.string "isomeric_smiles"
    t.string "inchi"
    t.string "inchi_key"
    t.string "iupac_name"
    t.string "title"
    t.string "x_log_p"
    t.string "exact_mass"
    t.string "monoisotopic_mass"
    t.string "tpsa"
    t.string "complexity"
    t.string "charge"
    t.string "h_bond_donor_count"
    t.string "h_bond_acceptor_count"
    t.string "rotatable_bond_count"
    t.string "heavy_atom_count"
    t.string "isotope_atom_count"
    t.string "atom_stereo_count"
    t.string "defined_atom_stereo_count"
    t.string "undefined_atom_stereo_count"
    t.string "bond_stereo_count"
    t.string "defined_bond_stereo_count"
    t.string "undefined_bond_stereo_count"
    t.string "covalent_unit_count"
    t.string "patent_count"
    t.string "patent_family_count"
    t.string "literature_count"
    t.string "cid"
    t.index ["protein_id"], name: "index_inhibitors_on_protein_id"
  end

  create_table "intermediates", force: :cascade do |t|
    t.integer "product_number"
    t.integer "step_number"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "products", force: :cascade do |t|
    t.string "hazards"
    t.float "melting_point_c"
    t.float "boiling_point_c"
    t.string "reported_reactions"
    t.integer "reported_yeild"
    t.integer "actual_yeild"
    t.float "obtained_mass"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "proteins", force: :cascade do |t|
    t.string "title"
    t.string "description"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.boolean "oncogenic"
    t.boolean "allostery"
    t.integer "size_kb"
    t.string "gene"
  end

  create_table "reactants", force: :cascade do |t|
    t.string "reactive_functional_group"
    t.float "mass_used"
    t.float "concentration"
    t.string "hazards"
    t.float "melting_point_c"
    t.float "boiling_point_c"
    t.float "mass_available"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "reaction_schemes", force: :cascade do |t|
    t.integer "number_of_reactions"
    t.integer "number_of_reatants"
    t.string "similar_reaction_schemes"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  create_table "reactions", force: :cascade do |t|
    t.integer "number_of_products"
    t.integer "number_of_reatants"
    t.string "solvent"
    t.string "catalysts"
    t.integer "step_number"
    t.string "type_of_reaction"
    t.float "temperature_c"
    t.string "similar_reaction_schemes"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
  end

  add_foreign_key "inhibitors", "proteins"
end
