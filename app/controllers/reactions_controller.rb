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
              
            reactants= ["NC(=S)N(C)C", "[Cl]CC(=O)c1ccccc1"]
            rxn ="[NH2:1][C:2](=[S:3])[NH0:4].[C:5](=[O:6])[C:7][Cl:8]>>[N:4][c:2]1[s:3][c:5][c:7][n:1]1"
          
            @generated_products = RdkitService.get_predicted_products(reactants, rxn)
          end
    
          def show
            @reaction = Reaction.find(params[:id])
            reactants= ["NC(=S)N(C)C", "[Cl]CC(=O)c1ccccc1"]
            rxn ="[NH2:1][C:2](=[S:3])[NH0:4].[C:5](=[O:6])[C:7][Cl:8]>>[N:4][c:2]1[s:3][c:5][c:7][n:1]1"
          
            # @generated_products = RdkitService.get_predicted_products(reactants, rxn)
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