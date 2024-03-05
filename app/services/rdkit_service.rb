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

    def self.get_reaction_from_mols(rxn_dict)
        begin
            response = HTTParty.post('http://localhost:1500/get_rxn_from_mols', body: { rxn_dict: rxn_dict }.to_json, headers: { 'Content-Type' => 'application/json' })
        rescue URI::InvalidURIError
            return nil
        end 
    
        if response.code == 200
            reaction_smarts = response['reaction_smiles']
            puts "Received reaction SMILES: #{reaction_smarts}"
            return reaction_smarts
        else
            puts "Error: #{response.code}, #{response.body}"
            return nil
        end
    end
    
  def self.get_mols_from_rxn(rxn_smarts)
    begin
      response = HTTParty.post('http://localhost:1500/get_mols_from_rxn', body: { rxn_smarts: rxn_smarts }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      rxn_dict = JSON.parse(response.body)
      puts "Received reaction dictionary: #{rxn_dict}"
      return rxn_dict
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.get_standardized_smarts(rxn_smarts)
    begin
        response = HTTParty.post('http://localhost:1500/get_standardized_smarts', body: { rxn_smarts: rxn_smarts }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
        return nil
    end 

    if response.code == 200
        standardized_smarts = response['standardized_smarts']
        puts "Received standardized SMARTS: #{standardized_smarts}"
        return standardized_smarts
    else
        puts "Error: #{response.code}, #{response.body}"
        return nil
    end
end

end

