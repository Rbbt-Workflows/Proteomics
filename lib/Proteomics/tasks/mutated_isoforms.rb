require 'rbbt/sources/organism'
require 'Proteomics/identifiers'
require 'Proteomics/neighbours'
require 'Proteomics/annotators'
require 'Proteomics/unfold'

module Proteomics

  helper :parse_mi do |mi, organism=Organism.default_code("Hsa")|
    isoform, residue = nil

    case
    when (m = mi.match(/^([^:]*):([A-Z])(\d+)([A-Z])(?:#{FOLD_SEP}.*)?$/))
      next if m[2] == m[4]
      isoform = m[1]
      residue = m[3].to_i
    when (m = mi.match(/^([^:]*):(\d+)(?:#{FOLD_SEP}.*)?$/))
      isoform = m[1]
      residue = m[2].to_i
    when mi.match(/:UTR/)
    when (m = mi.match(/^([^:]*):([A-Z\*])(\d+)([A-Z\*]+)(?:#{FOLD_SEP}.*)?$/i))
    when (m = mi.match(/^(.*):([A-Z\*])(\d+)(?:FrameShift|Indel)\(([+\-A-Z*]*)\)(?:#{FOLD_SEP}.*)?$/))
      isoform = m[1]
      residue = m[3].to_i
      substitution = m[4]
      #if substitution[-1] == "*"
      #  residue = [residue, -1] * ":"
      #else
      #  residue = [residue.to_s, (residue + substitution.length - 1).to_s] * ":"
      #end
    else
      raise "Unknown mutated isoform: #{mi}"
    end

    isoform = Proteomics.gene2isoform(isoform, organism) if isoform && isoform !~ /^ENS[A-Z]*P/

    [isoform, residue]
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :only_pdb, :boolean, "Only consider PDB neighbours", false
  input :just_one, :boolean, "Consider only neighbours from first PDB that contains them", true
  input :use_i3d, :boolean, "Use also i3d models", true
  task :mi_neighbours => :tsv do |mis,organism,only_pdb,just_one,use_i3d|

    cpus =  config :cpus, :mi_neighbours, :neighbours, :proteomics, :Proteomics, :default => 2

    annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "PDB", "Neighbours"], :type => :double, :namespace => organism
    annotations.init
    TSV.traverse mis, :cpus => cpus, :bar => self.progress_bar("Mutated Isoform neighbours"), :into => annotations, :type => :array do |mi|

      isoform, residue = parse_mi mi, organism

      next if isoform.nil?

      n = Proteomics.neighbours(isoform, [residue], organism, 5, only_pdb, just_one, use_i3d)

      next if n.empty?

      pdbs = []
      ns = []
      n.each do |r,v|
        Misc.zip_fields(v).each do |_p,_pos,pdb,n|
          next if n.empty?
          pdbs << pdb
          ns << n.split(";")
        end
      end

      next if ns.empty?

      [mi, [residue, pdbs.compact.uniq, ns.flatten.compact.uniq]]
    end
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :distance, :float, "Distance with partner residue", 8
  task :mi_interfaces => :tsv do |mis,organism,distance|
    cpus =  config :cpus, :mi_interfaces, :interfaces, :proteomics, :Proteomics, :default => 2

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "Partner (Ensembl Protein ID)", "PDB", "Partner Residues"], :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse mis, :cpus => cpus, :bar => self.progress_bar("Mutated Isoform interfaces"), :into => mi_annotations, :type => :array do |mi|

      isoform, residue = parse_mi mi, organism

      next if isoform.nil?

      n = Proteomics.interface_neighbours_i3d(isoform.dup, [residue], organism, distance)

      next if n.nil? or n.empty?

      all_annots = []
      n.each do |r,v|
        _pos,part,pdb,n = v
        next if part.nil? or part.empty?
        all_annots << [residue, part, pdb, n]
      end

      [mi, Misc.zip_fields(all_annots)]
    end
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true, :required => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  task :annotate_mi => :tsv do |mis,organism,database|

    annotator = ANNOTATORS[database]
    raise ParameterException, "Database not identified: #{ database }" if annotator.nil?
    annotator.organism = organism

    cpus = config :cpus, :annotate_mi, :annotate, :proteomics, :Proteomics, :default => 2

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => annotator.fields, :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse mis, :cpus => cpus, :bar => self.progress_bar("Annot. #{ database }"), :into => mi_annotations, :type => :array do |mi|

      isoform, residue = parse_mi(mi)
      
      next if isoform.nil?

      annotations = annotator.annotate isoform, residue, organism
      next if annotations.nil?

      [mi, annotations]
    end
  end



  dep :mi_neighbours
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  input :simplify, :boolean, "Simplify output lists", false
  task :annotate_mi_neighbours => :tsv do |organism,database,simplify|

    stream = Misc.open_pipe do |sin|
      dumper = TSV::Dumper.new :key_fields => "Mutated Isoform", :fields =>  ["Neighbours"],:type => :flat, :namespace => organism
      dumper.init
      TSV.traverse step(:mi_neighbours), :into => dumper do |mi,values|
        mi = mi.first if Array === mi
        isoform = mi.split(":").first
        residue, pdbs, ns = values
        neighbours = ns.collect{|pos| [isoform, pos] * ":"}
        neighbours.extend MultipleResult
        [mi, neighbours]
      end
      Misc.consume_stream(dumper.stream, false, sin)
      sin.close
    end

    Proteomics.unfold_traverse(stream, Proteomics, :annotate_mi, :mutated_isoforms, :database => database, :organism => organism, :unfold_field => "Neighbour", :simplify => simplify) do |mi,nmi,values|
      values[0] = nmi.split(":").last.scan(/\d+/).first
      values.collect{|l| Array === l ? l * "|" : l }
    end
  end


end
