require File.join(File.expand_path(File.dirname(__FILE__)), '../..', 'test_helper.rb')
require 'rbbt/PDB/isoform'

class TestIsoform < Test::Unit::TestCase
  def test_uniprot
    assert PDB.isoform_uniprot_pdbs("ENSP00000257189").any?
  end
  
  def test_i3d_models
    Log.severity = 0
    ppp PDB.isoform_i3d_pdbs("ENSP00000257189")
    assert PDB.isoform_i3d_pdbs("ENSP00000257189").any?
    assert PDB.isoform_i3d_pdbs("ENSP00000257189").keys.first.include?("5eqx")
  end

  def test_i3d_interactions
    ppp PDB.isoform_i3d_interaction_pdb("ENSP00000252725")
    assert PDB.isoform_i3d_interaction_pdb("ENSP00000252725").select("Partner (Ensembl Protein ID)" => "ENSP00000380431").any?
    assert PDB.isoform_i3d_interaction_pdb("ENSP00000352918").select("Partner (Ensembl Protein ID)" => "ENSP00000252725").any?
  end
end

