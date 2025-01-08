require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestInteractome3d < Test::Unit::TestCase
  def _test_pdbs_proteins
    pdbs = Interactome3d.proteins.tsv
    uni = pdbs.keys[1000]
    pdb = pdbs[uni]["FILENAME"].first

    iii pdb
    ppp Open.read(Interactome3d.pdb(pdb))
  end

  def test_pdbs_interactions
    pdbs = Interactome3d.interactions.tsv
    uni = pdbs.keys[1000]
    pdb = pdbs[uni]["FILENAME"].first

    iii pdb
    ppp Open.read(Interactome3d.pdb(pdb))
  end
end

