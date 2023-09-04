class ReactionsController < ApplicationController
        def index
            @reactions = Reactions.all
        end

        def new
        end
    
        def create
            reaction_scheme = Reaction.new({
                number_of_reactions: params[:reaction_scheme][:number_of_reactions],
                number_of_reatants: params[:reaction_scheme][:number_of_reatants],
                similar_reaction_schemes: params[:reaction_scheme][:similar_reaction_schemes]
              })
          
            reaction.save
          
            redirect_to '/reactions'
          end
    
          def show
            @reaction = Reaction.find(params[:id])
        
          end
    
          def edit
            @reaction = Reaction.find(params[:id])
          end
    
          def update
            @reaction = Reaction.find(params[:id])
            reaction.update({
                number_of_reactions: params[:reaction_scheme][:number_of_reactions],
                number_of_reatants: params[:reaction_scheme][:number_of_reatants],
                similar_reaction_schemes: params[:reaction_scheme][:similar_reaction_schemes]
              })
              reaction.save
            redirect_to "/reactions/#{reaction.id}"
          end
    
          def destroy
            Reaction.destroy(params[:id])
            redirect_to '/reactions'
          end
          
      end