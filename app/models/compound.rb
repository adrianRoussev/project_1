class Compound < ApplicationRecord
    belongs_to :compoundable, polymorphic: true, optional: true,  dependent: :destroy
  
    has_many :reactants, dependent: :destroy
    has_many :products, dependent: :destroy
    has_many :catalysts, dependent: :destroy
  
    has_many :reactions_from_reactants, through: :reactants, source: :reaction, dependent: :destroy
    has_many :reactions_from_products, through: :products, source: :reaction, dependent: :destroy
    has_many :reactions_from_catalysts, through: :catalysts, source: :reaction, dependent: :destroy




    def get_params_from_rdkit(parsed_data, compound)

        parsed_data = RdkitService.get_rdkit_params(compound.smiles)
        parsed_data_hash = parsed_data.to_h

        parsed_data_hash.each do |key, value|
            param = key.downcase.to_sym
            compound.update(param => parsed_data_hash[key])
            compound.save
        end
    end
                #     avgipc: parsed_data['AvgIpc'],
                #     bcut2d_chghi: parsed_data['BCUT2D_CHGHI'],
                #     bcut2d_chglo: parsed_data['BCUT2D_CHGLO'],
                #     bcut2d_logphi: parsed_data['BCUT2D_LOGPHI'],
                #     bcut2d_logplow: parsed_data['BCUT2D_LOGPLOW'],
                #     bcut2d_mrhi: parsed_data['BCUT2D_MRHI'],
                #     bcut2d_mrlo: parsed_data['BCUT2D_MRLOW'],
                #     bcut2d_mwhi: parsed_data['BCUT2D_MWHI'],
                #     bcut2d_mwlow: parsed_data['BCUT2D_MWLOW'],
                #     balabanj: parsed_data['BalabanJ'],
                #     bertzct: parsed_data['BertzCT'],
                #     chi0: parsed_data['Chi0'],
                #     chi0n: parsed_data['Chi0n'],
                #     chi0v: parsed_data['Chi0v'],
                #     chi1: parsed_data['Chi1'],
                #     chi1n: parsed_data['Chi1n'],
                #     chi1v: parsed_data['Chi1v'],
                #     chi2n: parsed_data['Chi2n'],
                #     chi2v: parsed_data['Chi2v'],
                #     chi3n: parsed_data['Chi3n'],
                #     chi3v: parsed_data['Chi3v'],
                #     chi4n: parsed_data['Chi4n'],
                #     chi4v: parsed_data['Chi4v'],
                #     estate_vsa1: parsed_data['EState_VSA1'],
                #     estate_vsa10: parsed_data['EState_VSA10'],
                #     estate_vsa11: parsed_data['EState_VSA11'],
                #     estate_vsa2: parsed_data['EState_VSA2'],
                #     estate_vsa3: parsed_data['EState_VSA3'],
                #     estate_vsa4: parsed_data['EState_VSA4'],
                #     estate_vsa5: parsed_data['EState_VSA5'],
                #     estate_vsa6: parsed_data['EState_VSA6'],
                #     estate_vsa7: parsed_data['EState_VSA7'],
                #     estate_vsa8: parsed_data['EState_VSA8'],
                #     estate_vsa9: parsed_data['EState_VSA9'],
                #     exactmolwt: parsed_data['ExactMolWt'],
                #     fpdensitymorgan1: parsed_data['FpDensityMorgan1'],
                #     fpdensitymorgan2: parsed_data['FpDensityMorgan2'],
                #     fpdensitymorgan3: parsed_data['FpDensityMorgan3'],
                #     fractioncsp3: parsed_data['FractionCSP3'],
                #     hallkieralpha: parsed_data['HallKierAlpha'],
                #     heavyatomcount: parsed_data['HeavyAtomCount'],
                #     heavyatommolwt: parsed_data['HeavyAtomMolWt'],
                #     ipc: parsed_data['Ipc'],
                #     kappa1: parsed_data['Kappa1'],
                #     kappa2: parsed_data['Kappa2'],
                #     kappa3: parsed_data['Kappa3'],
                #     labuteasa: parsed_data['LabuteASA'],
                #     maxabsestateindex: parsed_data['MaxAbsEStateIndex'],
                #     maxabspartialcharge: parsed_data['MaxAbsPartialCharge'],
                #     maxestateindex: parsed_data['MaxEStateIndex'],
                #     maxpartialcharge: parsed_data['MaxPartialCharge'],
                #     minabsestateindex: parsed_data['MinAbsEStateIndex'],
                #     minabspartialcharge: parsed_data['MinAbsPartialCharge'],
                #     minestateindex: parsed_data['MinEStateIndex'],
                #     minpartialcharge: parsed_data['MinPartialCharge'],
                #     mollogp: parsed_data['MolLogP'],
                #     molmr: parsed_data['MolMR'],
                #     molwt: parsed_data['MolWt'],
                #     nhohcount: parsed_data['NHOHCount'],
                #     nocount: parsed_data['NOCount'],
                #     numaliphaticcarbocycles: parsed_data['NumAliphaticCarbocycles'],
                #     numaliphaticheterocycles: parsed_data['NumAliphaticHeterocycles'],
                #     numaliphaticrings: parsed_data['NumAliphaticRings'],
                #     numaromaticcarbocycles: parsed_data['NumAromaticCarbocycles'],
                #     numaromaticheterocycles: parsed_data['NumAromaticHeterocycles'],
                #     numaromaticrings: parsed_data['NumAromaticRings'],
                #     numhacceptors: parsed_data['NumHAcceptors'],
                #     numhdonors: parsed_data['NumHDonors'],
                #     numheteroatoms: parsed_data['NumHeteroatoms'],
                #     numradicalelectrons: parsed_data['NumRadicalElectrons'],
                #     numrotatablebonds: parsed_data['NumRotatableBonds'],
                #     numsaturatedcarbocycles: parsed_data['NumSaturatedCarbocycles'],
                #     numsaturatedheterocycles: parsed_data['NumSaturatedHeterocycles'],
                #     numsaturatedrings: parsed_data['NumSaturatedRings'],
                #     numvalenceelectrons: parsed_data['NumValenceElectrons'],
                #     peoe_vsa1: parsed_data['PEOE_VSA1'],
                #     peoe_vsa10: parsed_data['PEOE_VSA10'],
                #     peoe_vsa11: parsed_data['PEOE_VSA11'],
                #     peoe_vsa12: parsed_data['PEOE_VSA12'],
                #     peoe_vsa13: parsed_data['PEOE_VSA13'],
                #     peoe_vsa14: parsed_data['PEOE_VSA14'],
                #     peoe_vsa2: parsed_data['PEOE_VSA2'],
                #     peoe_vsa3: parsed_data['PEOE_VSA3'],
                #     peoe_vsa4: parsed_data['PEOE_VSA4'],
                #     peoe_vsa5: parsed_data['PEOE_VSA5'],
                #     peoe_vsa6: parsed_data['PEOE_VSA6'],
                #     peoe_vsa7: parsed_data['PEOE_VSA7'],
                #     peoe_vsa8: parsed_data['PEOE_VSA8'],
                #     peoe_vsa9: parsed_data['PEOE_VSA9'],
                #     ringcount: parsed_data['RingCount'],
                #     smr_vsa1: parsed_data['SMR_VSA1'],
                #     smr_vsa10: parsed_data['SMR_VSA10'],
                #     smr_vsa2: parsed_data['SMR_VSA2'],
                #     smr_vsa3: parsed_data['SMR_VSA3'],
                #     smr_vsa4: parsed_data['SMR_VSA4'],
                #     smr_vsa5: parsed_data['SMR_VSA5'],
                #     smr_vsa6: parsed_data['SMR_VSA6'],
                #     smr_vsa7: parsed_data['SMR_VSA7'],
                #     smr_vsa8: parsed_data['SMR_VSA8'],
                #     smr_vsa9: parsed_data['SMR_VSA9'],
                #     slogp_vsa1: parsed_data['SlogP_VSA1'],
                #     slogp_vsa10: parsed_data['SlogP_VSA10'],
                #     slogp_vsa11: parsed_data['SlogP_VSA11'],
                #     slogp_vsa12: parsed_data['SlogP_VSA12'],
                #     slogp_vsa2: parsed_data['SlogP_VSA2'],
                #     slogp_vsa3: parsed_data['SlogP_VSA3'],
                #     slogp_vsa4: parsed_data['SlogP_VSA4'],
                #     slogp_vsa5: parsed_data['SlogP_VSA5'],
                #     slogp_vsa6: parsed_data['SlogP_VSA6'],
                #     slogp_vsa7: parsed_data['SlogP_VSA7'],
                #     slogp_vsa8: parsed_data['SlogP_VSA8'],
                #     slogp_vsa9: parsed_data['SlogP_VSA9'],
                #     tpsa: parsed_data['TPSA'],
                #     vsa_estate1: parsed_data['VSA_EState1'],
                #     vsa_estate10: parsed_data['VSA_EState10'],
                #     vsa_estate2: parsed_data['VSA_EState2'],
                #     vsa_estate3: parsed_data['VSA_EState3'],
                #     vsa_estate4: parsed_data['VSA_EState4'],
                #     vsa_estate5: parsed_data['VSA_EState5'],
                #     vsa_estate6: parsed_data['VSA_EState6'],
                #     vsa_estate7: parsed_data['VSA_EState7'],
                #     vsa_estate8: parsed_data['VSA_EState8'],
                #     vsa_estate9: parsed_data['VSA_EState9'],
                #     fr_al_coo: parsed_data['fr_Al_COO'],
                #     fr_al_oh: parsed_data['fr_Al_OH'],
                #     fr_al_oh_notert: parsed_data['fr_Al_OH_noTert'],
                #     fr_arn: parsed_data['fr_ArN'],
                #     fr_ar_coo: parsed_data['fr_Ar_COO'],
                #     fr_ar_n: parsed_data['fr_Ar_N'],
                #     fr_ar_nh: parsed_data['fr_Ar_NH'],
                #     fr_ar_oh: parsed_data['fr_Ar_OH'],
                #     fr_coo: parsed_data['fr_COO'],
                #     fr_coo2: parsed_data['fr_COO2'],
                #     fr_c_o: parsed_data['fr_C_O'],
                #     fr_c_o_nocoo: parsed_data['fr_C_O_noCOO'],
                #     fr_c_s: parsed_data['fr_C_S'],
                #     fr_hoccn: parsed_data['fr_HOCCN'],
                #     fr_imine: parsed_data['fr_Imine'],
                #     fr_nho: parsed_data['fr_NH0'],
                #     fr_nh1: parsed_data['fr_NH1'],
                #     fr_nh2: parsed_data['fr_NH2'],
                #     fr_n_o: parsed_data['fr_N_O'],
                #     fr_ndealkylation1: parsed_data['fr_Ndealkylation1'],
                #     fr_ndealkylation2: parsed_data['fr_Ndealkylation2'],
                #     fr_nhpyrrole: parsed_data['fr_Nhpyrrole'],
                #     fr_sh: parsed_data['fr_SH'],
                #     fr_aldehyde: parsed_data['fr_aldehyde'],
                #     fr_alkyl_carbamate: parsed_data['fr_alkyl_carbamate'],
                #     fr_alkyl_halide: parsed_data['fr_alkyl_halide'],
                #     fr_allylic_oxid: parsed_data['fr_allylic_oxid'],
                #     fr_amide: parsed_data['fr_amide'],
                #     fr_amidine: parsed_data['fr_amidine'],
                #     fr_aniline: parsed_data['fr_aniline'],
                #     fr_aryl_methyl: parsed_data['fr_aryl_methyl'],
                #     fr_azide: parsed_data['fr_azide'],
                #     fr_azo: parsed_data['fr_azo'],
                #     fr_barbitur: parsed_data['fr_barbitur'],
                #     fr_benzene: parsed_data['fr_benzene'],
                #     fr_benzodiazepine: parsed_data['fr_benzodiazepine'],
                #     fr_bicyclic: parsed_data['fr_bicyclic'],
                #     fr_diazo: parsed_data['fr_diazo'],
                #     fr_dihydropyridine: parsed_data['fr_dihydropyridine'],
                #     fr_epoxide: parsed_data['fr_epoxide'],
                #     fr_ester: parsed_data['fr_ester'],
                #     fr_ether: parsed_data['fr_ether'],
                #     fr_furan: parsed_data['fr_furan'],
                #     fr_guanido: parsed_data['fr_guanido'],
                #     fr_halogen: parsed_data['fr_halogen'],
                #     fr_hdrzine: parsed_data['fr_hdrzine'],
                #     fr_hdrzone: parsed_data['fr_hdrzone'],
                #     fr_imidazole: parsed_data['fr_imidazole'],
                #     fr_imide: parsed_data['fr_imide'],
                #     fr_isocyan: parsed_data['fr_isocyan'],
                #     fr_isothiocyan: parsed_data['fr_isothiocyan'],
                #     fr_ketone: parsed_data['fr_ketone'],
                #     fr_ketone_topliss: parsed_data['fr_ketone_Topliss'],
                #     fr_lactam: parsed_data['fr_lactam'],
                #     fr_lactone: parsed_data['fr_lactone'],
                #     fr_methoxy: parsed_data['fr_methoxy'],
                #     fr_morpholine: parsed_data['fr_morpholine'],
                #     fr_nitrile: parsed_data['fr_nitrile'],
                #     fr_nitro: parsed_data['fr_nitro'],
                #     fr_nitro_arom: parsed_data['fr_nitro_arom'],
                #     fr_nitro_arom_nonortho: parsed_data['fr_nitro_arom_nonortho'],
                #     fr_nitroso: parsed_data['fr_nitroso'],
                #     fr_oxazole: parsed_data['fr_oxazole'],
                #     fr_oxime: parsed_data['fr_oxime'],
                #     fr_para_hydroxylation: parsed_data['fr_para_hydroxylation'],
                #     fr_phenol: parsed_data['fr_phenol'],
                #     fr_phenol_noOrthoHbond: parsed_data['fr_phenol_noOrthoHbond'],
                #     fr_phos_acid: parsed_data['fr_phos_acid'],
                #     fr_phos_ester: parsed_data['fr_phos_ester'],
                #     fr_piperdine: parsed_data['fr_piperdine'],
                #     fr_piperzine: parsed_data['fr_piperzine'],
                #     fr_priamide: parsed_data['fr_priamide'],
                #     fr_prisulfonamd: parsed_data['fr_prisulfonamd'],
                #     fr_pyridine: parsed_data['fr_pyridine'],
                #     fr_quatN: parsed_data['fr_quatN'],
                #     fr_sulfide: parsed_data['fr_sulfide'],
                #     fr_sulfonamd: parsed_data['fr_sulfonamd'],
                #     fr_sulfone: parsed_data['fr_sulfone'],
                #     fr_term_acetylene: parsed_data['fr_term_acetylene'],
                #     fr_tetrazole: parsed_data['fr_tetrazole'],
                #     fr_thiazole: parsed_data['fr_thiazole'],
                #     fr_thiocyan: parsed_data['fr_thiocyan'],
                #     fr_thiophene: parsed_data['fr_thiophene'],
                #     fr_unbrch_alkane: parsed_data['fr_unbrch_alkane'],
                #     fr_urea: parsed_data['fr_urea'],
                #     qed: parsed_data['qed']
                # })

    #         compound.save
          
      
    # end
end


# def update_attributes_from_parsed_data(parsed_data)
#     parsed_data.each do |parsed_key, value|
#     model_key = convert_to_model_key(parsed_key)

#     if self.respond_to?("#{model_key}=")
#         self.send("#{model_key}=", value)
#     end
  
#   end
# end

# private

# def convert_to_model_key(parsed_key)
#   parsed_key.downcase.to_sym
# end

# end


#cont
#   my_model_instance = MyModel.new
# parsed_data = { 'Key1' => 'value1', 'Key2' => 'value2' }

# my_model_instance.update_attributes_from_parsed_data(parsed_data)
# my_model_instance.save


# #   compound = Compound.create(molecular_weight: 100, molecular_formula: 'C6H12O6')

# # # Creating an inhibitor associated with a compound
# # inhibitor = Inhibitor.create(name: 'Some Inhibitor')
# # inhibitor.compoundable = compound
# # inhibitor.save

# # # Accessing compounds associated with an inhibitor
# # inhibitor.compoundable

# # # Accessing inhibitors associated with a compound
# # compound.inhibitors


