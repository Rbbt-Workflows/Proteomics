require File.join(File.expand_path(File.dirname(__FILE__)), '../..', 'test_helper.rb')
require 'rbbt/PDB/distances'
require 'Proteomics/neighbours'

class TestDistances < Test::Unit::TestCase
  def test_interaction_surfaces
    tsv = PDB.isoform_i3d_interaction_pdb("ENSP00000252725")
    url = tsv[tsv.keys.first]["URL"].first
    p = PDB.interface_residues(url).first
    assert_not_equal p.first.split(":").first, p.last.split(":").first
  end
end

