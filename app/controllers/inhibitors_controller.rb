class InhibitorsController < ApplicationController
    def index
      @inhibitors = Inhibitor.all
    end
  
    def show
      @inhibitor = Inhibitor.find(params[:id])
      smiles = @inhibitor.smiles
      name = @inhibitor.name
      @properties = PubChemService.import_compound_properties(smiles, name)
    
      if @properties.present?
    @inhibitor.create_or_update_from_csv_data(@properties,@inhibitor)

      end
      # @inhibitor = Inhibitor.find_by(smiles: smiles)
      # render json: @inhibitor, serializer: InhibitorSerializer
    end
  
    def edit
        @inhibitor = Inhibitor.find(params[:id])
    end
    
    def update
        inhibitor = Inhibitor.find(params[:id])
        if inhibitor.smiles == nil 
      inhibitor.update({
        structure: params[:inhibitor][:structure],
        name: params[:inhibitor][:name],
        target_site: params[:inhibitor][:target_site],
        smiles: params[:inhibitor][:smiles]
          })
      inhibitor.save
      else 
        inhibitor.update({
          structure: params[:inhibitor][:structure],
          name: params[:inhibitor][:name],
          target_site: params[:inhibitor][:target_site],
          smiles: params[:inhibitor][:smiles]
            })
        inhibitor.save
        smiles = inhibitor.smiles
        name = inhibitor.name


        @properties = PubChemService.import_compound_properties(smiles, name)
       
        inhibitor.create_or_update_from_csv_data(@properties, inhibitor) 
     
        inhibitor.save

        end
    
      python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
    image_file_path = Rails.root.join('public', "#{inhibitor.name}.png")

    system("python3 #{python_script_path} '#{inhibitor.smiles}' '#{image_file_path}'")

    inhibitor.update(image_path: "/#{inhibitor.name}.png")



      inhibitor.save
  
      redirect_to "/inhibitors/#{inhibitor.id}"
    end

    def new_compoundable
      @inhibitor = Inhibitor.find(params[:id])
    end

    # def create_compoundable
    #   @inhibitor = Inhibitor.find(params[:id])
    #   @inhibitor_compoundable = Compound.create({:compoundable => @inhibitor,
    #       name: params[:compound][:name],
    #       molecular_weight: params[:compound][:molecular_weight],
    #       molecular_formula: params[:compound][:molecular_formula],
    #       smiles: params[:compound][:smiles],
    #       nucleophiles: params[:compound][:nucleophiles],
    #       electrophiles: params[:compound][:electrophiles],
    #       charge: params[:compound][:charge],
    #       log_p: params[:compound][:log_p],
    #       mass_available: params[:compound][:mass_available]
    #       })
          # compoundable_id: params[:compound][:compoundable][:compoundable_id],
        # compoundable_type: params[:compound][:compoundable][:compoundable_type]
        # .compoundable({compoundable_id: params[:compoundable][:compoundable_id],
        #   compoundable_type: params[:compoundable][:compoundable_type]})
      
    # end
    

    
    def create_compoundable
      @inhibitor = Inhibitor.find(params[:id])
      @inhibitor_compoundable = @inhibitor.compounds.create({
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
      
        redirect_to "/compounds"
    end
    
    
    def destroy
      Inhibitor.destroy(params[:id])
      
    end
  end


  