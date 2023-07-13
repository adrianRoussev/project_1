class ProteinsController < ApplicationController
    def index
      @proteins = Protein.all
    end

    def new
    end

    def create
        protein = Protein.new({
          title: params[:protein][:title],
          description: params[:protein][:description]
          })
      
        protein.save
      
        redirect_to '/proteins'
      end

      def show
        @protein = Protein.find(params[:id])
      end

      def edit
        @protein = Protein.find(params[:id])
      end

      def update
        protein = Protein.find(params[:id])
        protein.update({
          title: params[:protein][:title],
          description: params[:protein][:description]
          })
        protein.save
        redirect_to "/proteins/#{protein.id}"
      end

      def destroy
        Protein.destroy(params[:id])
        redirect_to '/proteins'
      end
      
  end