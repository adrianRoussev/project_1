class CompoundsController < ApplicationController
    def index
        @compounds = Compound.all
    end

    def new
    end

    def create
        compound = Compound.new({ 
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
           compound.update(compound_params)
                              
            redirect_to "/compounds"
        end
       

      def destroy
        Compound.destroy(params[:id])
        redirect_to '/compounds'
      end

      private

  def compound_params
    params.permit(:name, :molecular_weight, :molecular_formula, :smiles, :nucleophiles, :electrophiles, :charge, :log_p, :mass_available)
  end
end
      
  