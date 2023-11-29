class CompoundsController < ApplicationController
    def index
        @compounds = Compound.all
    end

    def new
    end
   
    def create
      inhibitor = Inhibitor.create
      compound = Compound.create({
        compoundable: inhibitor, 
        name: params[:compound][:name],
        molecular_weight: params[:compound][:molecular_weight],
        molecular_formula: params[:compound][:molecular_formula],
        smiles: params[:compound][:smiles],
        nucleophiles: params[:compound][:nucleophiles],
        electrophiles: params[:compound][:electrophiles],
        charge: params[:compound][:charge],
        log_p: params[:compound][:log_p],
        mass_available: params[:compound][:mass_available]
      })
              compound.image_path = "/#{compound.name}.png"
        compound.save
      
        redirect_to '/compounds'
      end

      def show
        @compound = Compound.find(params[:id])
      end

      def edit
        @compound = Compound.find(params[:id])
      end

      def update
          compound = Compound.find(params[:id])
          compound.update({ 
            name: params[:compound][:name],
            molecular_weight: params[:compound][:molecular_weight],
            molecular_formula: params[:compound][:molecular_formula],
            smiles: params[:compound][:smiles],
            nucleophiles: params[:compound][:nucleophiles],
            electrophiles: params[:compound][:electrophiles],
            charge: params[:compound][:charge],
            log_p: params[:compound][:log_p],
            mass_available: params[:compound][:mass_available]
          })
          compound.image_path = "/#{compound.name}.png"
          compound.save
            redirect_to "/compounds"

            python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
            image_file_path = Rails.root.join('public', "#{compound.name}.png")
        
            system("python3 #{python_script_path} '#{compound.smiles}' '#{image_file_path}'")
        
            compound.update(image_path: "/#{compound.name}.png")
        
            compound.save

        end


        def new_inhibitor
          @proteins = Protein.all
          @compound = Compound.find(params[:id])
    
        end

        def create_inhibitor
          compound = Compound.find(params[:id])
          protein = Protein.find(params[:protein_id])
          inhibitor = Inhibitor.create({
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
            protein_id: params[:protein_id]
          })
inhibitor.save

compound.update(compoundable: inhibitor)

          
        end
  
        #     smiles = compound.smiles
        #     name = compound.name
  
        #   @properties = PubChemService.import_compound_properties(smiles, name)
       
        #   inhibitor.create_or_update_from_csv_data(@properties, inhibitor) 
     
        # inhibitor.save
         
     

      def destroy
        Compound.destroy(params[:id])
        redirect_to '/compounds'
      end

      # private
      # def compound_params
      #   params.permit(
      #     :name,
      #     :molecular_weight,
      #     :molecular_formula,
      #     :smiles,
      #     :nucleophiles,
      #     :electrophiles,
      #     :charge,
      #     :log_p,
      #     :mass_available
      #   )
      # end
      

end
      
  