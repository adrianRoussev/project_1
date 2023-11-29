class CompoundInhibitorsController < ApplicationController
    def new
      @compound = Compound.find(params[:id])

    end


    # def create
    #   compound = Compound.find(params[:compound_id])
    #   inhibitor = Inhibitor.create({
    #     name: params[:inhibitor][:name],
    #     structure: params[:inhibitor][:structure],
    #     lipinskis: params[:inhibitor][:lipinskis],
    #     target_site: params[:inhibitor][:target_site],
    #     no_h_donnors: params[:inhibitor][:no_h_donnors],
    #     no_h_acceptors: params[:inhibitor][:no_h_acceptors],
    #     no_rotatable_bonds: params[:inhibitor][:no_rotatable_bonds],
    #     no_heavy_atoms: params[:inhibitor][:no_heavy_atoms],
    #     smiles: params[:inhibitor][:smiles],
    #     molecular_formula: params[:inhibitor][:molecular_formula],
    #     molecular_weight: params[:inhibitor][:molecular_weight],
    #     canonical_smiles: params[:inhibitor][:canonical_smiles],
    #     isomeric_smiles: params[:inhibitor][:isomeric_smiles],
    #     inchi: params[:inhibitor][:inchi],
    #     inchi_key: params[:inhibitor][:inchi_key],
    #     iupac_name: params[:inhibitor][:iupac_name],
    #     title: params[:inhibitor][:title],
    #     x_log_p: params[:inhibitor][:x_log_p],
    #     exact_mass: params[:inhibitor][:exact_mass],
    #     monoisotopic_mass: params[:inhibitor][:monoisotopic_mass],
    #     tpsa: params[:inhibitor][:tpsa],
    #     complexity: params[:inhibitor][:complexity],
    #     charge: params[:inhibitor][:charge],
    #     h_bond_donor_count: params[:inhibitor][:h_bond_donor_count],
    #     h_bond_acceptor_count: params[:inhibitor][:h_bond_acceptor_count],
    #     rotatable_bond_count: params[:inhibitor][:rotatable_bond_count],
    #     heavy_atom_count: params[:inhibitor][:heavy_atom_count],
    #     isotope_atom_count: params[:inhibitor][:isotope_atom_count],
    #     atom_stereo_count: params[:inhibitor][:atom_stereo_count],
    #     defined_atom_stereo_count: params[:inhibitor][:defined_atom_stereo_count],
    #     undefined_atom_stereo_count: params[:inhibitor][:undefined_atom_stereo_count],
    #     bond_stereo_count: params[:inhibitor][:bond_stereo_count],
    #     defined_bond_stereo_count: params[:inhibitor][:defined_bond_stereo_count],
    #     undefined_bond_stereo_count: params[:inhibitor][:undefined_bond_stereo_count],
    #     covalent_unit_count: params[:inhibitor][:covalent_unit_count],
    #     patent_count: params[:inhibitor][:patent_count],
    #     patent_family_count: params[:inhibitor][:patent_family_count],
    #     literature_count: params[:inhibitor][:literature_count],
    #     compound_id: params[:id]
    #   })
    
    
  
    #   python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
    #   image_file_path = Rails.root.join('public', "#{inhibitor.name}.png")
  
    #   system("python3 #{python_script_path} '#{inhibitor.smiles}' '#{image_file_path}'")
  
    #   inhibitor.update(image_path: "/#{inhibitor.name}.png")
  
  
    #    inhibitor.save
    #     redirect_to "/inhibitors"
     
    # end

  
    def create
        @inhibitor = Inhibitor.new({
          name: params[:inhibitor][:name],
          structure: params[:inhibitor][:structure],
          lipinskis: params[:inhibitor][:lipinskis],
          target_site: params[:inhibitor][:target_site],
          no_h_donnors: params[:inhibitor][:no_h_donnors],
          no_h_acceptors: params[:inhibitor][:no_h_acceptors],
          no_rotatable_bonds: params[:inhibitor][:no_rotatable_bonds],
          no_heavy_atoms: params[:inhibitor][:no_heavy_atoms],
          smiles: params[:inhibitor][:smiles],
          molecular_formula: params[:inhibitor][:molecular_formula],
          molecular_weight: params[:inhibitor][:molecular_weight],
          canonical_smiles: params[:inhibitor][:canonical_smiles],
          isomeric_smiles: params[:inhibitor][:isomeric_smiles],
          inchi: params[:inhibitor][:inchi],
          inchi_key: params[:inhibitor][:inchi_key],
          iupac_name: params[:inhibitor][:iupac_name],
          title: params[:inhibitor][:title],
          x_log_p: params[:inhibitor][:x_log_p],
          exact_mass: params[:inhibitor][:exact_mass],
          monoisotopic_mass: params[:inhibitor][:monoisotopic_mass],
          tpsa: params[:inhibitor][:tpsa],
          complexity: params[:inhibitor][:complexity],
          charge: params[:inhibitor][:charge],
          h_bond_donor_count: params[:inhibitor][:h_bond_donor_count],
          h_bond_acceptor_count: params[:inhibitor][:h_bond_acceptor_count],
          rotatable_bond_count: params[:inhibitor][:rotatable_bond_count],
          heavy_atom_count: params[:inhibitor][:heavy_atom_count],
          isotope_atom_count: params[:inhibitor][:isotope_atom_count],
          atom_stereo_count: params[:inhibitor][:atom_stereo_count],
          defined_atom_stereo_count: params[:inhibitor][:defined_atom_stereo_count],
          undefined_atom_stereo_count: params[:inhibitor][:undefined_atom_stereo_count],
          bond_stereo_count: params[:inhibitor][:bond_stereo_count],
          defined_bond_stereo_count: params[:inhibitor][:defined_bond_stereo_count],
          undefined_bond_stereo_count: params[:inhibitor][:undefined_bond_stereo_count],
          covalent_unit_count: params[:inhibitor][:covalent_unit_count],
          patent_count: params[:inhibitor][:patent_count],
          patent_family_count: params[:inhibitor][:patent_family_count],
          literature_count: params[:inhibitor][:literature_count],
        })
        # inhibitor.compoundable = compound

        # inhibitor.name = compound.name
        # inhibitor.structure = compound.molecular_formula
        # inhibitor.smiles = compound.smiles
        # inhibitor.image_path = compound.image_path
        
       

          if @inhibitor.save
            redirect_to "/inhibitors", notice: 'Inhibitor created successfully.'
          end

          smiles = @inhibitor.smiles
          name = @inhibitor.name

        @properties = PubChemService.import_compound_properties(smiles, name)
     
        @inhibitor.create_or_update_from_csv_data(@properties, @inhibitor) 
   
      @inhibitor.save
       
      end
    end

     #   compound = Compound.create(molecular_weight: 100, molecular_formula: 'C6H12O6')

# # Creating an inhibitor associated with a compound
# inhibitor = Inhibitor.create(name: 'Some Inhibitor')
# inhibitor.compoundable = compound
# inhibitor.save

# # Accessing compounds associated with an inhibitor
# inhibitor.compoundable

# # Accessing inhibitors associated with a compound
# compound.inhibitors