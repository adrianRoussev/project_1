# app/services/pub_chem_service.rb
class PubChemService
    include HTTParty
    base_uri 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
  
    def self.import_compound_properties(smiles, name)
      # Compose the URL to request the first 29 properties for the given CID
      url = "/compound/smiles/#{smiles}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES," \
            "InChI,InChIKey,IUPACName,Title,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge," \
            "HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount," \
            "AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount," \
            "DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,PatentCount," \
            "PatentFamilyCount,LiteratureCount/CSV"

    begin response = get(url)
    
    rescue URI::InvalidURIError
            nil
    return response

    end 

    if response.success?
    
        csv_data_raw = response.body

        if csv_data_raw.present?
       
            csv_data_parsed = CSV.parse(csv_data_raw, headers: true)
            return csv_data_parsed
     
        else
            return nil

        end

    else
        return nil
    end
  end
end
