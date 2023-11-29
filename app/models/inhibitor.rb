
    # require 'pycall/import'
    class Inhibitor < ApplicationRecord
      belongs_to :protein
      include Compoundable
      has_many :reactions
  
   

      def create_or_update_from_csv_data(csv_data_parsed, inhibitor)
       if csv_data_parsed == nil 
       nil
       else 
        csv_data_parsed.each do |row|
          inhibitor.assign_attributes({
            cid: row['CID'],
            molecular_formula: row['MolecularFormula'],
            molecular_weight: row['MolecularWeight'],
            canonical_smiles: row['CanonicalSMILES'],
            isomeric_smiles: row['IsomericSMILES'],
            inchi: row['InChI'],
            inchi_key: row['InChIKey'],
            iupac_name: row['IUPACName'],
            title: row['Title'],
            x_log_p: row['XLogP'],
            exact_mass: row['ExactMass'],
            monoisotopic_mass: row['MonoisotopicMass'],
            tpsa: row['TPSA'],
            complexity: row['Complexity'],
            charge: row['Charge'],
            h_bond_donor_count: row['HBondDonorCount'],
            h_bond_acceptor_count: row['HBondAcceptorCount'],
            rotatable_bond_count: row['RotatableBondCount'],
            heavy_atom_count: row['HeavyAtomCount'],
            isotope_atom_count: row['IsotopeAtomCount'],
            atom_stereo_count: row['AtomStereoCount'],
            defined_atom_stereo_count: row['DefinedAtomStereoCount'],
            undefined_atom_stereo_count: row['UndefinedAtomStereoCount'],
            bond_stereo_count: row['BondStereoCount'],
            defined_bond_stereo_count: row['DefinedBondStereoCount'],
            undefined_bond_stereo_count: row['UndefinedBondStereoCount'],
            covalent_unit_count: row['CovalentUnitCount'],
            patent_count: row['PatentCount'],
            patent_family_count: row['PatentFamilyCount'],
            literature_count: row['LiteratureCount']
          })
    
          inhibitor.save
          end
        end
      end
    end
    