class ReactionSchemesController < ApplicationController
        def index
            @reaction_scheme = ReactionScheme.all
        end

        def new
        end
    
        def create
            reaction_scheme = ReactionScheme.new({
                number_of_reactions: params[:reaction_scheme][:number_of_reactions],
                number_of_reatants: params[:reaction_scheme][:number_of_reatants],
                similar_reaction_schemes: params[:reaction_scheme][:similar_reaction_schemes]
              })
          
            reaction_scheme.save
          
            redirect_to '/proteins'
          end
    
          def show
            @reaction_scheme = ReactionScheme.find(params[:id])
        
          end
    
          def edit
            @reaction_scheme = ReactionScheme.find(params[:id])
          end
    
          def update
            @reaction_scheme = ReactionScheme.find(params[:id])
            reaction_scheme.update({
                number_of_reactions: params[:reaction_scheme][:number_of_reactions],
                number_of_reatants: params[:reaction_scheme][:number_of_reatants],
                similar_reaction_schemes: params[:reaction_scheme][:similar_reaction_schemes]
              })
              reaction_scheme.save
            redirect_to "/proteins/#{protein.id}"
          end
    
          def destroy
            ReactionScheme.destroy(params[:id])
            redirect_to '/proteins'
          end
          
      end