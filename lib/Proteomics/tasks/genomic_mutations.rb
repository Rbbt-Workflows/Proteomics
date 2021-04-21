module Proteomics

  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  input :principal, :boolean, "Consider only principal isoforms", true
  input :simplify, :boolean, "Simplify output lists", false
  dep Sequence, :mutated_isoforms_fast, :principal => :principal, :non_synonymous => true
  task :annotate_dna => :tsv do |organism,database,principal,simplify|
    Proteomics.unfold_traverse(step(:mutated_isoforms_fast), Proteomics, :annotate_mi, :mutated_isoforms, :database => database, :organism => organism, :unfold_field => "Mutated Isoform", :key_field => "Genomic Mutation", :simplify => simplify) do |mut,mi,values|
      values.collect{|l| Array === l ? l * "|" : l }
    end
  end

  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :principal, :boolean, "Consider only principal isoforms", true
  dep Sequence, :mutated_isoforms_fast, :principal => :principal, :non_synonymous => true
  task :dna_neighbours => :tsv do |organism,principal|
    Proteomics.unfold_traverse(step(:mutated_isoforms_fast), Proteomics, :mi_neighbours, :mutated_isoforms, :organism => organism, :unfold_field => "Mutated Isoform", :key_field => "Genomic Mutation") do |mut,mi,values|
      values.collect{|l| Array === l ? l * "|" : l }
    end
  end

  dep :dna_neighbours
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  input :simplify, :boolean, "Simplify output lists", false
  task :annotate_dna_neighbours => :tsv do |organism,database,simplify|

    stream = Misc.open_pipe do |sin|
      dumper = TSV::Dumper.new :key_fields => "Mutated Isoform", :fields =>  ["Neighbours"],:type => :flat, :namespace => organism
      dumper.init
      TSV.traverse step(:dna_neighbours), :into => dumper do |mut,values|
        mut = mut.first if Array === mut
        res = Misc.zip_fields(values).collect do |mi,res,pdb,ns|
          isoform = mi.split(":").first
          neighbours = ns.split(";").collect{|pos| [isoform, pos] * ":"}
          [mut, neighbours]
        end

        res.extend MultipleResult

        res
      end
      Misc.consume_stream(dumper.stream, false, sin)
      sin.close
    end

    Proteomics.unfold_traverse(stream, Proteomics, :annotate_mi, :mutated_isoforms, :database => database, :organism => organism, :unfold_field => "Neighbour", :key_field => "Genomic Mutation", :simplify => simplify) do |mi,nmi,values|
      values[0] = nmi.split(":").last.scan(/\d+/).first
      values.collect{|l| Array === l ? l * "|" : l }
    end
  end




  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  input :principal, :boolean, "Consider only principal isoforms", true
  input :neighbours, :boolean, "Annotate neighbours", false
  input :simplify, :boolean, "Simplify output lists", false
  dep Sequence, :mutated_isoforms_fast, :principal => :principal, :non_synonymous => true
  task :dna_interfaces => :tsv do |organism,database,principal,neighbours,simplify|
    Proteomics.unfold_traverse(step(:mutated_isoforms_fast), Proteomics, :mi_interfaces, :mutated_isoforms, :database => database, :organism => organism, :unfold_field => "Mutated Isoform", :key_field => "Genomic Mutation", :simplify => simplify) do |mut,mi,values|
      values.collect{|l| Array === l ? l * "|" : l }
    end
  end
end
