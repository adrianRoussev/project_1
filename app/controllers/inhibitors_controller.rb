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
        structure: params["Chemical Structure"],
        name: params["Chemical Name"],
        target_site: params["Target Site"],
        smiles: params["Smiles"]
          })
      inhibitor.save
      else 
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

    def destroy
      Inhibitor.destroy(params[:id])
      redirect_to "/proteins"
    end
  end


  