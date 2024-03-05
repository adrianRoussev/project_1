class ReactionsController < ApplicationController
        def index
            @all_reactions = Reaction.all
            @inhibitor = Inhibitor.find(params[:id])
            @reactions = Reaction.where(inhibitor_id: params[:id])

        end

        def index_all
          inhibitor_id = params[:inhibitor_id]
          @reactions = Reaction.where(inhibitor_id: inhibitor_id)
          @inhibitors = Inhibitor.all
        end

        def search
          inhibitor_id = params[:inhibitor_id]
          @reactions = Reaction.where(inhibitor_id: inhibitor_id)
          @inhibitors = Inhibitor.all
        end


        def new
          @compounds= Compound.all
          @inhibitor = Inhibitor.find(params[:inhibitor_id])
          @inhibitor_reactions=@inhibitor.reactions.all
        end
    
        def create
          @compounds= Compound.all
          @inhibitor = Inhibitor.find(params[:inhibitor_id])
          @reaction = @inhibitor.reactions.build({
            solvent: params[:reaction][:solvent],
            step_number: params[:reaction][:step_number],
            reaction_smiles: params[:reaction][:reaction_smiles],
            temperature_c: params[:reaction][:temperature_c],
            inhibitor_id: params[:id],
            reactants_attributes: map_compounds(params[:reaction][:reactants]) ,
            products_attributes: map_compounds(params[:reaction][:products]),
            catalysts_attributes: map_compounds(params[:reaction][:catalysts])
  })
            


            reactants = @reaction.reactants
            products = @reaction.products
            catalysts = @reaction.catalysts

            reaction_smiles = @reaction.create_reaction(products, reactants, catalysts)

            reaction_smarts = RdkitService.get_standardized_smarts(reaction_smiles)

            @reaction.update(reaction_smiles: reaction_smarts)
    
            flash[:success] = 'Reaction created successfully.'

            # Make below a method in a helper file #

            ###############
            python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
            image_file_path = Rails.root.join('public', "#{@inhibitor.id}-#{@reaction.id}.png")
        
            system("python3 #{python_script_path} '#{reaction_smarts}' '#{image_file_path}'")
        
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
            @compounds= Compound.all
            # @inhibitor = Inhibitor.find(params[:inhibitor_id])
            @reaction = Reaction.find(params[:id])
            @reaction_id = @reaction.id
            
            reactants = @reaction.reactants
            products = @reaction.products
            catalysts = @reaction.catalysts
          end
    
    
          def update
            @reaction = Reaction.find(params[:id])
            @inhibitor = Inhibitor.find(@reaction.inhibitor_id)

            if @reaction.update(reaction_params)

              #update reactions
              @reaction.reactants.destroy_all
              reactants_array = params[:reaction][:reactants].map { |r| { compound_id: r } }
              @reaction.update(reactants_attributes: reactants_array)
        
              # Update products
              @reaction.products.destroy_all
              products_array = params[:reaction][:products].map { |r| { compound_id: r } }
              @reaction.update(products_attributes: products_array)
  
              # Update catalysts
              @reaction.catalysts.destroy_all
              catalysts_array = params[:reaction][:catalysts].map { |r| { compound_id: r } }
              @reaction.update(catalysts_attributes: catalysts_array)
  

              reactants = @reaction.reactants
              products = @reaction.products
              catalysts = @reaction.catalysts
  
              reaction_smiles = @reaction.create_reaction(products, reactants, catalysts)
              reaction_smarts = RdkitService.get_standardized_smarts(reaction_smiles)
            
              @reaction.update(reaction_smiles: reaction_smarts)
              flash[:success] = 'Reaction created successfully.'
  
              generate_image_for_reaction(@inhibitor, @reaction) 
              flash[:success] = 'Reaction updated successfully.'
      
              redirect_to "/inhibitors/#{@inhibitor.id}/reactions"
            
            else
              
              flash[:error] = 'Failed to update reaction.'
              render :edit
            
            end
          
          end
  
          def destroy
            Reaction.destroy(params[:id])
            # redirect_to '/reactions'
          end

            
          private
  
        def reaction_params
          params.require(:reaction).permit(:solvent, :step_number, :temperature_c, reactants: [:id, :compound_id], products: [:id, :compound_id], catalysts: [:id, :compound_id], reactants_attributes: [:compound_id], products_attributes: [:compound_id], catalysts_attributes: [:compound_id])
        end


        def generate_image_for_reaction(inhibitor, reaction)
          
          python_script_path = Rails.root.join('python_scripts', 'smiles_render.py')
          image_file_path = Rails.root.join('public', "#{inhibitor.id}-#{reaction.id}.png")
          reaction_smiles = reaction.reaction_smiles
          
          system("python3 #{python_script_path} '#{reaction_smiles}' '#{image_file_path}'")
        
        end


        def map_compounds(compounds)
          return [] if compounds.blank?
        
          compounds.map { |c| { compound_id: c } }
        end
          
      end