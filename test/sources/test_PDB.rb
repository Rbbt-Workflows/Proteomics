require File.join(File.expand_path(File.dirname(__FILE__)), '..', 'test_helper.rb')
require 'sources/PDB'

class TestPDBHelper < Test::Unit::TestCase
  def test_stream
    stream = PDB.pdb_stream('3dzy')
    assert stream.read.include? "ATOM"
  end

  def test_chains
    chains = PDB.pdb_chain_sequences('3dzy')
    assert_equal "H", chains["A"][132]
    assert_equal "I", chains["A"][133]
  end
end

