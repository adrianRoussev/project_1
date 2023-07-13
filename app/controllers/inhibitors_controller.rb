class InhibitorsController < ApplicationController
    def index
      @inhibitors = Inhibitor.all
    end
  
    def show
      @inhibitor = Inhibitor.find(params[:id])
    end
  
    def edit
        @inhibitor = Inhibitor.find(params[:id])
    end
    
    def update
        inhibitor = Inhibitor.find(params[:id])
      inhibitor.update({
        structure: params["Chemical Structure"],
        name: params["Chemical Name"],
        lipinski: params["Lipinski violation"],
        target_site: params["active-site or allosteric"],
        weight: params["molecular weight (Da)"]
      })
      inhibitor.save
  
      redirect_to "/inhibitors/#{inhibitors.id}"
    end
  end