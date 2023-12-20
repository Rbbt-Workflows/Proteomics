require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestAlphafold < Test::Unit::TestCase
  def test_json
    uniprot = "P00520"
    info = AlphaFold.json(uniprot)
    assert_equal "Abl1", info[0]['gene']
  end

  def test_pdbs
    uniprot = "P00520"
    list = AlphaFold.pdbs(uniprot)
    iii list
  end
end

