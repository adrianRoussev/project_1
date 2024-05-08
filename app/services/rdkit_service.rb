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

  def self.get_compounds_from_smarts(rxn_smarts, compound_type)
    begin
      response = HTTParty.post('http://localhost:1500/get_compounds_from_smarts', body: { rxn_smarts: rxn_smarts, compound_type: compound_type }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    return response if response.nil?

    if response.code == 200
      compounds = JSON.parse(response.body)['compounds']
      puts "Received compounds: #{compounds}"
      return compounds
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.get_mol_functional_groups(smiles_string)
    begin
      response = HTTParty.post('http://localhost:1500/get_mol_functional_groups', body: { smiles_string: smiles_string }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      functional_groups = JSON.parse(response.body)['functional_groups']
      puts "Received functional groups: #{functional_groups}"
      return functional_groups
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.find_reaction_functional_groups(smarts)
    begin
      response = HTTParty.post('http://localhost:1500/find_reaction_functional_groups', body: { smarts: smarts }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      functional_groups = JSON.parse(response.body)['functional_groups']
      puts "Received functional groups: #{functional_groups}"
      return functional_groups
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.find_changing_functional_groups(smarts)
    begin
      response = HTTParty.post('http://localhost:1500/find_changing_functional_groups', body: { smarts: smarts }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      changing_groups = JSON.parse(response.body)['changing_groups']
      puts "Received changing functional groups: #{changing_groups}"
      return changing_groups
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.get_base_rxn(smarts, reagents = false)
    begin
      response = HTTParty.post('http://localhost:1500/get_base_rxn', body: { smarts: smarts, reagents: reagents }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      base_rxn_smarts = response['base_rxn_smarts']
      puts "Received base reaction SMILES: #{base_rxn_smarts}"
      return base_rxn_smarts
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end

  def self.run_base_rxn(smarts, reagents = false)
    begin
      response = HTTParty.post('http://localhost:1500/run_base_rxn', body: { smarts: smarts, reagents: reagents }.to_json, headers: { 'Content-Type' => 'application/json' })
    rescue URI::InvalidURIError
      return nil
    end 

    if response.code == 200
      base_rxn_results = JSON.parse(response.body)['base_rxn_results']
      puts "Received base reaction results: #{base_rxn_results}"
      return base_rxn_results
    else
      puts "Error: #{response.code}, #{response.body}"
      return nil
    end
  end
end

