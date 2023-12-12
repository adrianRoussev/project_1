require 'json'

class CompoundsController < ApplicationController
    def index
        @compounds = Compound.all

    end

    def new
    end
   
    def create
      # inhibitor = Inhibitor.create
      compound = Compound.create({
        # compoundable: inhibitor, 
        name: params[:compound][:name],
        molecular_formula: params[:compound][:molecular_formula],
        smiles: params[:compound][:smiles],
        nucleophiles: params[:compound][:nucleophiles],
        electrophiles: params[:compound][:electrophiles],
        charge: params[:compound][:charge],
        log_p: params[:compound][:log_p],
        mass_available: params[:compound][:mass_available],
        molwt: params[:compound][:molwt]

      })
      compound.image_path = "/#{compound.name}.png"
      compound.get_params_from_rdkit(@molecular_params, compound)

      compound.save

      python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
      image_file_path = Rails.root.join('public', "#{compound.name}.png")
        
      system("python3 #{python_script_path} '#{compound.smiles}' '#{image_file_path}'")
            
      compound.update(image_path: "/#{compound.name}.png")
      
      compound.save
            
      redirect_to '/compounds'
      end

      def show
        @compound = Compound.find(params[:id])
        
        @compound.get_params_from_rdkit(@molecular_params, @compound)

        @compound.save

      end

      def edit
        @compound = Compound.find(params[:id])
        @molecular_params = RdkitService.get_rdkit_params(@compound.smiles)

        @value = @compound.get_value_from_rdkit(@molecular_params)
     

      end

      def update
          @compound = Compound.find(params[:id])
          inhibitor = Inhibitor.create
           @compound.update({ 
            compoundable: inhibitor,
            name: params[:compound][:name],
            molecular_formula: params[:compound][:molecular_formula],
            smiles: params[:compound][:smiles],
            nucleophiles: params[:compound][:nucleophiles],
            electrophiles: params[:compound][:electrophiles],
            charge: params[:compound][:charge],
            log_p: params[:compound][:log_p],
            mass_available: params[:compound][:mass_available],
            molwt: params[:compound][:molwt]
            
          })
          @compound.save

         
          
          

          @compound.image_path = "/#{@compound.name}.png"
          @compound.save
            

            python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
            image_file_path = Rails.root.join('public', "#{@compound.name}.png")
        
            system("python3 #{python_script_path} '#{@compound.smiles}' '#{image_file_path}'")
            
            @compound.update(image_path: "/#{@compound.name}.png")
            @compound.save


            redirect_to '/compounds', notice: 'Compound was successfully updated.'

        
           
        end


        def new_inhibitor
          @proteins = Protein.all
          @compound = Compound.find(params[:id])
    
        end

        #   create an inhibitor then set the compound to belong to THAT SPECIFIC inhibitor
        #   note:- a polymorphic compound belonging to an inhibitor can change its type and 
        #          become the child of a different type 
        #        - a certain type that belongs to a compound cannot
        #   Currently in a reaction we are creating types belonging to the compound and where
        #   types cannot be changed
        #  
        #   *questions*  1. Should we be creating the types in the reverse? ie. compound belonging to the type? 
        #                2. should both be an option? 
        #   currently for Reaction, we create types using a compound ID ie. Parent: Compound, Child: Type, Type has a compound ID
        #   
        #   could we instead just assign the compondable type and create 

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
    # inhibitor.compoundable.update(compoundable_type: 'Inhibitor', compoundable_id: @inhibitor.id)

          
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
      
  