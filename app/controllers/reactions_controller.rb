class ReactionsController < ApplicationController
        def index
            @all_reactions = Reaction.all
            @inhibitor = Inhibitor.find(params[:id])
            @reactions = Reaction.where(inhibitor_id: params[:id])
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
            solvent: params[:reaction][:solvent],
            step_number: params[:reaction][:step_number],
            reaction_smiles: params[:reaction][:reaction_smiles],
            temperature_c: params[:reaction][:temperature_c],
            inhibitor_id: params[:id],
            reactants_attributes: params[:reaction][:reactants].map { |r| { compound_id: r } },
            products_attributes: params[:reaction][:products].map { |s| { compound_id: s  } },
            catalysts_attributes: params[:reaction][:catalysts].map { |c| { compound_id: c } }
  })
            @reaction.save
            
           


                  
            reactants = @reaction.reactants
            products = @reaction.products
            catalysts = @reaction.catalysts
            reaction_smiles = @reaction.create_reaction(products, reactants, catalysts)

            @reaction.update(reaction_smiles: reaction_smiles)

            @reaction.save

            # Make below a method in a helper file #

            ###############
            python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
            image_file_path = Rails.root.join('public', "#{@inhibitor.id}-#{@reaction.id}.png")
        
            system("python3 #{python_script_path} '#{reaction_smiles}' '#{image_file_path}'")
        
            @reaction.update(image_path: "/#{@inhibitor.id}-#{@reaction.id}.png")
        
            @reaction.save

            redirect_to "/inhibitors/#{@inhibitor.id}/reactions"
            ################

            
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