require 'Proteomics/identifiers'
require 'Proteomics/neighbours'
require 'rbbt/sources/organism'

module Proteomics

  helper :parse_mi do |mi|
    isoform, residue = nil

    case
    when (m = mi.match(/^(.*):([A-Z])(\d+)([A-Z])$/))
      next if m[2] == m[4]
      isoform = m[1]
      residue = m[3].to_i
    when (m = mi.match(/^(.*):(\d+)$/))
      isoform = m[1]
      residue = m[2].to_i
    when mi.match(/:UTR/)
    when (m = mi.match(/^(.*):([A-Z\*])(\d+)([A-Z\*]+)$/i))
    else
      raise "Unknown mutated isoform: #{mi}"
    end

    if isoform[0..3] != "ENSP"
      mi = Proteomics.name2pi(mi, organism)
      next if mi.nil?
      isoform = mi.split(":").first
    end if isoform

    [isoform, residue]
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :mi_neighbours => :tsv do |mis,organism|

    cpus =  config :cpus, :mi_neighbours, :proteomics, :neighbours, :Proteomics, :default => 2

    annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "PDB", "Neighbours"], :type => :double, :namespace => organism
    annotations.init
    TSV.traverse mis, :cpus => cpus, :bar => self.progress_bar("Mutated Isoform neighbours"), :into => annotations, :type => :array do |mi|

      isoform, residue = parse_mi mi

      next if isoform.nil?

      n = Proteomics.neighbours(isoform, [residue], organism)

      next if n.empty?

      pdbs = []
      ns = []
      n.each do |r,v|
        _p,_pos,pdb,n = v
        next if n.empty?
        pdbs << pdb
        ns << n.split(";")
      end

      next if ns.empty?

      [mi, [residue, pdbs, ns]]
    end
  end
end
