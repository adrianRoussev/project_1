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

    # def self.get_raw_rdkit_data(smiles)

    #     begin response = HTTParty.post('http://localhost:5000/get_molecule_params', body: { smiles: smiles }.to_json, headers: { 'Content-Type' => 'application/json' })
    
    #     rescue URI::InvalidURIError
    #             nil
        
    #     return response

    #     end 

    #     if response.code == 200

    #         rdkit_data_raw = response.body
    #         rdkit_data_string= rdkit_data_string
    #         puts "Received molecule parameters: #{rdkit_data_raw}"
    #         return rdkit_data_string

    #     else

    #         puts "Error: #{response.code}, #{response.body}"
        
    #     end
    # end


end

