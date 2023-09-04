# app/serializers/inhibitor_serializer.rb
class InhibitorObject
  def Initialize 
    
      
      MolecularFormula: inhibitor.molecular_formula,
      MolecularWeight: inhibitor.molecular_weight,
      CanonicalSMILES: inhibitor.canonical_smiles,
      IsomericSMILES: inhibitor.isomeric_smiles,
      InChI: inhibitor.inchi,
      InChIKey: inhibitor.inchi_key,
      IUPACName: inhibitor.iupac_name,
      Title: inhibitor.title,
      XLogP: inhibitor.x_log_p,
      ExactMass: inhibitor.exact_mass,
      MonoisotopicMass: inhibitor.monoisotopic_mass,
      TPSA: inhibitor.tpsa,
      Complexity: inhibitor.complexity,
      Charge: inhibitor.charge,
      HBondDonorCount: inhibitor.h_bond_donor_count,
      HBondAcceptorCount: inhibitor.h_bond_acceptor_count,
      RotatableBondCount: inhibitor.rotatable_bond_count,
      HeavyAtomCount: inhibitor.heavy_atom_count,
      IsotopeAtomCount: inhibitor.isotope_atom_count,
      AtomStereoCount: inhibitor.atom_stereo_count,
      DefinedAtomStereoCount: inhibitor.defined_atom_stereo_count,
      UndefinedAtomStereoCount: inhibitor.undefined_atom_stereo_count,
      BondStereoCount: inhibitor.bond_stereo_count,
      DefinedBondStereoCount: inhibitor.defined_bond_stereo_count,
      UndefinedBondStereoCount: inhibitor.undefined_bond_stereo_count,
      CovalentUnitCount: inhibitor.covalent_unit_count,
      PatentCount: inhibitor.patent_count,
      PatentFamilyCount: inhibitor.patent_family_count,
      LiteratureCount: inhibitor.literature_count
    }
  end
end
