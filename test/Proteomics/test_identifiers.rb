require File.join(File.expand_path(File.dirname(__FILE__)), '..', 'test_helper.rb')
require 'rbbt-util'
require 'Proteomics/identifiers'

class TestIdentifiers < Test::Unit::TestCase
  def test_ensp2ensp
    isoform = "ENSP00000337194" 

    assert_equal isoform, Proteomics.gene2isoform(isoform)
    assert_equal isoform, Proteomics.gene2isoform("PRPF4B")
    assert_equal nil, Proteomics.gene2isoform("NOVALIDNAME")
  end

  def test_name2pi
    name = "PRPF4B:F854C" 
    mi = "ENSP00000337194:F854C" 

    assert_equal mi, Proteomics.name2pi(name)
  end
end

