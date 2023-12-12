class RdkitService
    include HTTParty
    
    def self.get_rdkit_params(smiles)

        begin response = HTTParty.post('http://localhost:5000/get_molecule_params', body: { smiles: smiles }.to_json, headers: { 'Content-Type' => 'application/json' })
    
        rescue URI::InvalidURIError
                nil
        
        return response

        end 

        if response.code == 200

            rdkit_data_raw = response.body
            puts "Received molecule parameters: #{rdkit_data_raw}"

            if rdkit_data_raw.present?
   
                molecule_params_parsed = JSON.parse(rdkit_data_raw, headers: true)
                puts "Parsed molecule parameters: #{molecule_params_parsed}"
                return molecule_params_parsed
        
            else
        
                return nil

            end

        else

            puts "Error: #{response.code}, #{response.body}"
        
        end
    end


    def self.get_predicted_products(reactants, rxn)

        begin response = HTTParty.post('http://localhost:6000/get_rxn_products', body: { reactants: reactants, rxn: rxn}.to_json, headers: { 'Content-Type' => 'application/json' })
    
        rescue URI::InvalidURIError
                nil
        
        return response

        end 

        if response.code == 200

            rdkit_products_raw = response.body
            puts "Received products: #{rdkit_products_raw}"

            if rdkit_products_raw.present?
   
                products_parsed = JSON.parse(rdkit_products_raw, headers: true)
                puts "Parsed predicted products: #{products_parsed}"
                return products_parsed
        
            else
        
                return nil

            end

        else

            puts "Error: #{response.code}, #{response.body}"
        
        end
    end

end

