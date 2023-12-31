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

ActiveRecord::Schema[7.0].define(version: 2023_12_11_034610) do
  # These are extensions that must be enabled in order to support this database
  enable_extension "plpgsql"

  create_table "catalysts", force: :cascade do |t|
    t.bigint "reaction_id", null: false
    t.bigint "compound_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.float "mass_used"
    t.index ["compound_id"], name: "index_catalysts_on_compound_id"
    t.index ["reaction_id"], name: "index_catalysts_on_reaction_id"
  end

  create_table "compoundables", force: :cascade do |t|
    t.string "compound_type", null: false
    t.bigint "compound_id", null: false
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "compoundable_type"
    t.bigint "compoundable_id"
    t.index ["compound_type", "compound_id"], name: "index_compoundables_on_compound"
    t.index ["compoundable_type", "compoundable_id"], name: "index_compoundables_on_compoundable"
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
    t.string "image_path"
    t.string "compound_type"
    t.string "compoundable_type"
    t.bigint "compoundable_id"
    t.string "molecular_params"
    t.float "avgipc"
    t.float "bcut2d_chghi"
    t.float "bcut2d_chglo"
    t.float "bcut2d_logphi"
    t.float "bcut2d_logplow"
    t.float "bcut2d_mrhi"
    t.float "bcut2d_mwhi"
    t.float "bcut2d_mwlow"
    t.float "balabanj"
    t.float "bertzct"
    t.float "chi0"
    t.float "chi0n"
    t.float "chi0v"
    t.float "chi1"
    t.float "chi1n"
    t.float "chi1v"
    t.float "chi2n"
    t.float "chi2v"
    t.float "chi3n"
    t.float "chi3v"
    t.float "chi4n"
    t.float "chi4v"
    t.float "estate_vsa1"
    t.float "estate_vsa10"
    t.float "estate_vsa11"
    t.float "estate_vsa2"
    t.float "estate_vsa3"
    t.float "estate_vsa4"
    t.float "estate_vsa5"
    t.float "estate_vsa6"
    t.float "estate_vsa7"
    t.float "estate_vsa8"
    t.float "estate_vsa9"
    t.float "exactmolwt"
    t.float "fpdensitymorgan1"
    t.float "fpdensitymorgan2"
    t.float "fpdensitymorgan3"
    t.float "fractioncsp3"
    t.float "hallkieralpha"
    t.integer "heavyatomcount"
    t.float "heavyatommolwt"
    t.float "ipc"
    t.float "kappa1"
    t.float "kappa2"
    t.float "kappa3"
    t.float "labuteasa"
    t.float "maxabsestateindex"
    t.float "maxabspartialcharge"
    t.float "maxestateindex"
    t.float "maxpartialcharge"
    t.float "minabsestateindex"
    t.float "minabspartialcharge"
    t.float "minestateindex"
    t.float "minpartialcharge"
    t.float "mollogp"
    t.float "molmr"
    t.float "molwt"
    t.integer "nhohcount"
    t.integer "nocount"
    t.integer "numaliphaticcarbocycles"
    t.integer "numaliphaticheterocycles"
    t.integer "numaliphaticrings"
    t.integer "numaromaticcarbocycles"
    t.integer "numaromaticheterocycles"
    t.integer "numaromaticrings"
    t.integer "numhacceptors"
    t.integer "numhdonors"
    t.integer "numheteroatoms"
    t.integer "numradicalelectrons"
    t.integer "numrotatablebonds"
    t.integer "numsaturatedcarbocycles"
    t.integer "numsaturatedheterocycles"
    t.integer "numsaturatedrings"
    t.integer "numvalenceelectrons"
    t.float "peoe_vsa1"
    t.float "peoe_vsa10"
    t.float "peoe_vsa11"
    t.float "peoe_vsa12"
    t.float "peoe_vsa13"
    t.float "peoe_vsa14"
    t.float "peoe_vsa2"
    t.float "peoe_vsa3"
    t.float "peoe_vsa4"
    t.float "peoe_vsa5"
    t.float "peoe_vsa6"
    t.float "peoe_vsa7"
    t.float "peoe_vsa8"
    t.float "peoe_vsa9"
    t.integer "ringcount"
    t.float "smr_vsa1"
    t.float "smr_vsa10"
    t.float "smr_vsa2"
    t.float "smr_vsa3"
    t.float "smr_vsa4"
    t.float "smr_vsa5"
    t.float "smr_vsa6"
    t.float "smr_vsa7"
    t.float "smr_vsa8"
    t.float "smr_vsa9"
    t.float "slogp_vsa1"
    t.float "slogp_vsa10"
    t.float "slogp_vsa11"
    t.float "slogp_vsa12"
    t.float "slogp_vsa2"
    t.float "slogp_vsa3"
    t.float "slogp_vsa4"
    t.float "slogp_vsa5"
    t.float "slogp_vsa6"
    t.float "slogp_vsa7"
    t.float "slogp_vsa8"
    t.float "slogp_vsa9"
    t.float "tpsa"
    t.float "vsa_estate1"
    t.float "vsa_estate10"
    t.float "vsa_estate2"
    t.float "vsa_estate3"
    t.float "vsa_estate4"
    t.float "vsa_estate5"
    t.float "vsa_estate6"
    t.float "vsa_estate7"
    t.float "vsa_estate8"
    t.float "vsa_estate9"
    t.integer "fr_al_coo"
    t.integer "fr_al_oh"
    t.integer "fr_al_oh_notert"
    t.integer "fr_arn"
    t.integer "fr_ar_coo"
    t.integer "fr_ar_n"
    t.integer "fr_ar_nh"
    t.integer "fr_ar_oh"
    t.integer "fr_coo"
    t.integer "fr_coo2"
    t.integer "fr_c_o"
    t.integer "fr_c_o_nocoo"
    t.integer "fr_c_s"
    t.integer "fr_hoccn"
    t.integer "fr_imine"
    t.integer "fr_nh1"
    t.integer "fr_nh2"
    t.integer "fr_n_o"
    t.integer "fr_ndealkylation1"
    t.integer "fr_ndealkylation2"
    t.integer "fr_nhpyrrole"
    t.integer "fr_sh"
    t.integer "fr_aldehyde"
    t.integer "fr_alkyl_carbamate"
    t.integer "fr_alkyl_halide"
    t.integer "fr_allylic_oxid"
    t.integer "fr_amide"
    t.integer "fr_amidine"
    t.integer "fr_aniline"
    t.integer "fr_aryl_methyl"
    t.integer "fr_azide"
    t.integer "fr_azo"
    t.integer "fr_barbitur"
    t.integer "fr_benzene"
    t.integer "fr_benzodiazepine"
    t.integer "fr_bicyclic"
    t.integer "fr_diazo"
    t.integer "fr_dihydropyridine"
    t.integer "fr_epoxide"
    t.integer "fr_ester"
    t.integer "fr_ether"
    t.integer "fr_furan"
    t.integer "fr_guanido"
    t.integer "fr_halogen"
    t.integer "fr_hdrzine"
    t.integer "fr_hdrzone"
    t.integer "fr_imidazole"
    t.integer "fr_imide"
    t.integer "fr_isocyan"
    t.integer "fr_isothiocyan"
    t.integer "fr_ketone"
    t.integer "fr_ketone_topliss"
    t.integer "fr_lactam"
    t.integer "fr_lactone"
    t.integer "fr_methoxy"
    t.integer "fr_morpholine"
    t.integer "fr_nitrile"
    t.integer "fr_nitro"
    t.integer "fr_nitro_arom"
    t.integer "fr_nitro_arom_nonortho"
    t.integer "fr_nitroso"
    t.integer "fr_oxazole"
    t.integer "fr_oxime"
    t.integer "fr_para_hydroxylation"
    t.integer "fr_phenol"
    t.integer "fr_phos_acid"
    t.integer "fr_phos_ester"
    t.integer "fr_piperdine"
    t.integer "fr_piperzine"
    t.integer "fr_priamide"
    t.integer "fr_prisulfonamd"
    t.integer "fr_pyridine"
    t.integer "fr_sulfide"
    t.integer "fr_sulfonamd"
    t.integer "fr_sulfone"
    t.integer "fr_term_acetylene"
    t.integer "fr_tetrazole"
    t.integer "fr_thiazole"
    t.integer "fr_thiocyan"
    t.integer "fr_thiophene"
    t.integer "fr_unbrch_alkane"
    t.integer "fr_urea"
    t.float "qed"
    t.integer "fr_nh0"
    t.integer "fr_quatn"
    t.integer "fr_phenol_noorthohbond"
    t.float "bcut2d_mrlow"
    t.index ["compoundable_type", "compoundable_id"], name: "index_compounds_on_compoundable"
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
    t.bigint "reaction_id", null: false
    t.bigint "compound_id", null: false
    t.index ["compound_id"], name: "index_products_on_compound_id"
    t.index ["reaction_id"], name: "index_products_on_reaction_id"
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
    t.bigint "reaction_id", null: false
    t.bigint "compound_id", null: false
    t.index ["compound_id"], name: "index_reactants_on_compound_id"
    t.index ["reaction_id"], name: "index_reactants_on_reaction_id"
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
    t.integer "step_number"
    t.string "type_of_reaction"
    t.float "temperature_c"
    t.datetime "created_at", null: false
    t.datetime "updated_at", null: false
    t.string "target_product_name"
    t.string "reaction_smiles"
    t.string "generated_products", default: [], array: true
    t.bigint "inhibitor_id", null: false
    t.string "image_path"
    t.index ["inhibitor_id"], name: "index_reactions_on_inhibitor_id"
  end

  add_foreign_key "catalysts", "compounds"
  add_foreign_key "catalysts", "reactions"
  add_foreign_key "inhibitors", "proteins"
  add_foreign_key "products", "compounds"
  add_foreign_key "products", "reactions"
  add_foreign_key "reactants", "compounds"
  add_foreign_key "reactants", "reactions"
  add_foreign_key "reactions", "inhibitors"
end
