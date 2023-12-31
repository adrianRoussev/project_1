class AddColumnsToCompounds < ActiveRecord::Migration[7.0]
  def change
    add_column :compounds, :AvgIpc, :float
    add_column :compounds, :BCUT2D_CHGHI, :float
    add_column :compounds, :BCUT2D_CHGLO, :float
    add_column :compounds, :BCUT2D_LOGPHI, :float
    add_column :compounds, :BCUT2D_LOGPLOW, :float
    add_column :compounds, :BCUT2D_MRHI, :float
    add_column :compounds, :BCUT2D_MRLOW, :float
    add_column :compounds, :BCUT2D_MWHI, :float
    add_column :compounds, :BCUT2D_MWLOW, :float
    add_column :compounds, :BalabanJ, :float
    add_column :compounds, :BertzCT, :float
    add_column :compounds, :Chi0, :float
    add_column :compounds, :Chi0n, :float
    add_column :compounds, :Chi0v, :float
    add_column :compounds, :Chi1, :float
    add_column :compounds, :Chi1n, :float
    add_column :compounds, :Chi1v, :float
    add_column :compounds, :Chi2n, :float
    add_column :compounds, :Chi2v, :float
    add_column :compounds, :Chi3n, :float
    add_column :compounds, :Chi3v, :float
    add_column :compounds, :Chi4n, :float
    add_column :compounds, :Chi4v, :float
    add_column :compounds, :EState_VSA1, :float
    add_column :compounds, :EState_VSA10, :float
    add_column :compounds, :EState_VSA11, :float
    add_column :compounds, :EState_VSA2, :float
    add_column :compounds, :EState_VSA3, :float
    add_column :compounds, :EState_VSA4, :float
    add_column :compounds, :EState_VSA5, :float
    add_column :compounds, :EState_VSA6, :float
    add_column :compounds, :EState_VSA7, :float
    add_column :compounds, :EState_VSA8, :float
    add_column :compounds, :EState_VSA9, :float
    add_column :compounds, :ExactMolWt, :float
    add_column :compounds, :FpDensityMorgan1, :float
    add_column :compounds, :FpDensityMorgan2, :float
    add_column :compounds, :FpDensityMorgan3, :float
    add_column :compounds, :FractionCSP3, :float
    add_column :compounds, :HallKierAlpha, :float
    add_column :compounds, :HeavyAtomCount, :integer
    add_column :compounds, :HeavyAtomMolWt, :float
    add_column :compounds, :Ipc, :float
    add_column :compounds, :Kappa1, :float
    add_column :compounds, :Kappa2, :float
    add_column :compounds, :Kappa3, :float
    add_column :compounds, :LabuteASA, :float
    add_column :compounds, :MaxAbsEStateIndex, :float
    add_column :compounds, :MaxAbsPartialCharge, :float
    add_column :compounds, :MaxEStateIndex, :float
    add_column :compounds, :MaxPartialCharge, :float
    add_column :compounds, :MinAbsEStateIndex, :float
    add_column :compounds, :MinAbsPartialCharge, :float
    add_column :compounds, :MinEStateIndex, :float
    add_column :compounds, :MinPartialCharge, :float
    add_column :compounds, :MolLogP, :float
    add_column :compounds, :MolMR, :float
    add_column :compounds, :MolWt, :float
    add_column :compounds, :NHOHCount, :integer
    add_column :compounds, :NOCount, :integer
    add_column :compounds, :NumAliphaticCarbocycles, :integer
    add_column :compounds, :NumAliphaticHeterocycles, :integer
    add_column :compounds, :NumAliphaticRings, :integer
    add_column :compounds, :NumAromaticCarbocycles, :integer
    add_column :compounds, :NumAromaticHeterocycles, :integer
    add_column :compounds, :NumAromaticRings, :integer
    add_column :compounds, :NumHAcceptors, :integer
    add_column :compounds, :NumHDonors, :integer
    add_column :compounds, :NumHeteroatoms, :integer
    add_column :compounds, :NumRadicalElectrons, :integer
    add_column :compounds, :NumRotatableBonds, :integer
    add_column :compounds, :NumSaturatedCarbocycles, :integer
    add_column :compounds, :NumSaturatedHeterocycles, :integer
    add_column :compounds, :NumSaturatedRings, :integer
    add_column :compounds, :NumValenceElectrons, :integer
    add_column :compounds, :PEOE_VSA1, :float
    add_column :compounds, :PEOE_VSA10, :float
    add_column :compounds, :PEOE_VSA11, :float
    add_column :compounds, :PEOE_VSA12, :float
    add_column :compounds, :PEOE_VSA13, :float
    add_column :compounds, :PEOE_VSA14, :float
    add_column :compounds, :PEOE_VSA2, :float
    add_column :compounds, :PEOE_VSA3, :float
    add_column :compounds, :PEOE_VSA4, :float
    add_column :compounds, :PEOE_VSA5, :float
    add_column :compounds, :PEOE_VSA6, :float
    add_column :compounds, :PEOE_VSA7, :float
    add_column :compounds, :PEOE_VSA8, :float
    add_column :compounds, :PEOE_VSA9, :float
    add_column :compounds, :RingCount, :integer
    add_column :compounds, :SMR_VSA1, :float
    add_column :compounds, :SMR_VSA10, :float
    add_column :compounds, :SMR_VSA2, :float
    add_column :compounds, :SMR_VSA3, :float
    add_column :compounds, :SMR_VSA4, :float
    add_column :compounds, :SMR_VSA5, :float
    add_column :compounds, :SMR_VSA6, :float
    add_column :compounds, :SMR_VSA7, :float
    add_column :compounds, :SMR_VSA8, :float
    add_column :compounds, :SMR_VSA9, :float
    add_column :compounds, :SlogP_VSA1, :float
    add_column :compounds, :SlogP_VSA10, :float
    add_column :compounds, :SlogP_VSA11, :float
    add_column :compounds, :SlogP_VSA12, :float
    add_column :compounds, :SlogP_VSA2, :float
    add_column :compounds, :SlogP_VSA3, :float
    add_column :compounds, :SlogP_VSA4, :float
    add_column :compounds, :SlogP_VSA5, :float
    add_column :compounds, :SlogP_VSA6, :float
    add_column :compounds, :SlogP_VSA7, :float
    add_column :compounds, :SlogP_VSA8, :float
    add_column :compounds, :SlogP_VSA9, :float
    add_column :compounds, :TPSA, :float
    add_column :compounds, :VSA_EState1, :float
    add_column :compounds, :VSA_EState10, :float
    add_column :compounds, :VSA_EState2, :float
    add_column :compounds, :VSA_EState3, :float
    add_column :compounds, :VSA_EState4, :float
    add_column :compounds, :VSA_EState5, :float
    add_column :compounds, :VSA_EState6, :float
    add_column :compounds, :VSA_EState7, :float
    add_column :compounds, :VSA_EState8, :float
    add_column :compounds, :VSA_EState9, :float
    add_column :compounds, :fr_Al_COO, :integer
    add_column :compounds, :fr_Al_OH, :integer
    add_column :compounds, :fr_Al_OH_noTert, :integer
    add_column :compounds, :fr_ArN, :integer
    add_column :compounds, :fr_Ar_COO, :integer
    add_column :compounds, :fr_Ar_N, :integer
    add_column :compounds, :fr_Ar_NH, :integer
    add_column :compounds, :fr_Ar_OH, :integer
    add_column :compounds, :fr_COO, :integer
    add_column :compounds, :fr_COO2, :integer
    add_column :compounds, :fr_C_O, :integer
    add_column :compounds, :fr_C_O_noCOO, :integer
    add_column :compounds, :fr_C_S, :integer
    add_column :compounds, :fr_HOCCN, :integer
    add_column :compounds, :fr_Imine, :integer
    add_column :compounds, :fr_NH0, :integer
    add_column :compounds, :fr_NH1, :integer
    add_column :compounds, :fr_NH2, :integer
    add_column :compounds, :fr_N_O, :integer
    add_column :compounds, :fr_Ndealkylation1, :integer
    add_column :compounds, :fr_Ndealkylation2, :integer
    add_column :compounds, :fr_Nhpyrrole, :integer
    add_column :compounds, :fr_SH, :integer
    add_column :compounds, :fr_aldehyde, :integer
    add_column :compounds, :fr_alkyl_carbamate, :integer
    add_column :compounds, :fr_alkyl_halide, :integer
    add_column :compounds, :fr_allylic_oxid, :integer
    add_column :compounds, :fr_amide, :integer
    add_column :compounds, :fr_amidine, :integer
    add_column :compounds, :fr_aniline, :integer
    add_column :compounds, :fr_aryl_methyl, :integer
    add_column :compounds, :fr_azide, :integer
    add_column :compounds, :fr_azo, :integer
    add_column :compounds, :fr_barbitur, :integer
    add_column :compounds, :fr_benzene, :integer
    add_column :compounds, :fr_benzodiazepine, :integer
    add_column :compounds, :fr_bicyclic, :integer
    add_column :compounds, :fr_diazo, :integer
    add_column :compounds, :fr_dihydropyridine, :integer
    add_column :compounds, :fr_epoxide, :integer
    add_column :compounds, :fr_ester, :integer
    add_column :compounds, :fr_ether, :integer
    add_column :compounds, :fr_furan, :integer
    add_column :compounds, :fr_guanido, :integer
    add_column :compounds, :fr_halogen, :integer
    add_column :compounds, :fr_hdrzine, :integer
    add_column :compounds, :fr_hdrzone, :integer
    add_column :compounds, :fr_imidazole, :integer
    add_column :compounds, :fr_imide, :integer
    add_column :compounds, :fr_isocyan, :integer
    add_column :compounds, :fr_isothiocyan, :integer
    add_column :compounds, :fr_ketone, :integer
    add_column :compounds, :fr_ketone_Topliss, :integer
    add_column :compounds, :fr_lactam, :integer
    add_column :compounds, :fr_lactone, :integer
    add_column :compounds, :fr_methoxy, :integer
    add_column :compounds, :fr_morpholine, :integer
    add_column :compounds, :fr_nitrile, :integer
    add_column :compounds, :fr_nitro, :integer
    add_column :compounds, :fr_nitro_arom, :integer
    add_column :compounds, :fr_nitro_arom_nonortho, :integer
    add_column :compounds, :fr_nitroso, :integer
    add_column :compounds, :fr_oxazole, :integer
    add_column :compounds, :fr_oxime, :integer
    add_column :compounds, :fr_para_hydroxylation, :integer
    add_column :compounds, :fr_phenol, :integer
    add_column :compounds, :fr_phenol_noOrthoHbond, :integer
    add_column :compounds, :fr_phos_acid, :integer
    add_column :compounds, :fr_phos_ester, :integer
    add_column :compounds, :fr_piperdine, :integer
    add_column :compounds, :fr_piperzine, :integer
    add_column :compounds, :fr_priamide, :integer
    add_column :compounds, :fr_prisulfonamd, :integer
    add_column :compounds, :fr_pyridine, :integer
    add_column :compounds, :fr_quatN, :integer
    add_column :compounds, :fr_sulfide, :integer
    add_column :compounds, :fr_sulfonamd, :integer
    add_column :compounds, :fr_sulfone, :integer
    add_column :compounds, :fr_term_acetylene, :integer
    add_column :compounds, :fr_tetrazole, :integer
    add_column :compounds, :fr_thiazole, :integer
    add_column :compounds, :fr_thiocyan, :integer
    add_column :compounds, :fr_thiophene, :integer
    add_column :compounds, :fr_unbrch_alkane, :integer
    add_column :compounds, :fr_urea, :integer
    add_column :compounds, :qed, :float
  end
end
