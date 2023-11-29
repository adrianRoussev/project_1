class ReactionsController < ApplicationController
        def index
            @reactions = Reaction.all
            @inhibitor = Inhibitor.find(params[:id])
        end

        def new
          @compounds= Compound.all
          @inhibitor = Inhibitor.find(params[:inhibitor_id])
@inhibitor_reactions=@inhibitor.reactions.all
        end
    
        def create
          @compounds= Compound.all
          @inhibitor = Inhibitor.find(params[:inhibitor_id])
          @reaction = @inhibitor.reactions.create({
           reactants: params[:reaction][:reactants],
           products: params[:reaction][:products],
           solvent: params[:reaction][:solvent],
           catalysts: params[:reaction][:catalysts],
           step_number: params[:reaction][:step_number],
           temperature_c: params[:reaction][:temperature_c],
           inhibitor_id: params[:id]
            })
            @reaction.save
            redirect_to "/inhibitors/#{@inhibitor.id}/reactions"
             
          end

          def predict_products
            @reaction = Reaction.find(params[:id])

          end
    
          def show
            @reaction = Reaction.find(params[:id])
        
          end
    
          def edit
            @reaction = Reaction.find(params[:id])
          end
    
          def update
            @reaction = Reaction.find(params[:id])
            @reaction.update({
             
             solvent: params[:reaction][:solvent],
             catalysts: params[:reaction][:catalysts],
             step_number: params[:reaction][:step_number],
             temperature_c: params[:reaction][:temperature_c]
              })
              @reaction.save
            redirect_to "/inhibitors/:id/reactions"
          end
    
          def destroy
            Reaction.destroy(params[:id])
            redirect_to '/reactions'
          end
          
      end