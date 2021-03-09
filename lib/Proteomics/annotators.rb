
module Proteomics

  ANNOTATORS = IndiferentHash.setup({})

  class Annotator

    def self.cache(code)
      Rbbt.var.Proteomics.ANNOTATORS[code.to_s]
    end

    def cache(code)
      Annotator.cache code
    end

    attr_accessor :fields, :organism
    def initialize(*fields, &block)
      @fields = fields.flatten
      @block = block
      class << self; self end.send(:define_method, :annotate, &block)
    end

    def annotate(*args)
      block.call(*args).collect do |l| 
        case l
        when nil
          nil
        when Array
          l.collect{|e| e.gsub(";", ',') }
        else
          l.gsub(";", ',')
        end
      end
    end
  end

  def self.add_annotator(name, *rest, &block)
    annotator = Annotator.new(*rest, &block)
    ANNOTATORS[name] = annotator
  end
end

require 'Proteomics/annotators/uniprot'
require 'Proteomics/annotators/appris'
require 'Proteomics/annotators/InterPro'
require 'Proteomics/annotators/COSMIC'


