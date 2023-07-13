class ProteinInhibitorsController < ApplicationController
  def index
    @protein = Protein.find(params[:id])
    @inhibitors = if params[:sort] == "alpha"
      Inhibitor.where(protein_id: params[:id]).order(name: :asc)
    else
      Inhibitor.where(protein_id: params[:id])
    end
end

  def new
    @protein = Protein.find(params[:id])
  end

  def create
    inhibitor = Inhibitor.new({
      name: params[:inhibitor][:name],
      structure: params[:inhibitor][:structure],
      protein_id: params[:id]
    })

     inhibitor.save
      redirect_to "/proteins/#{inhibitor.protein_id}/inhibitors"
   
  end
end

    
      