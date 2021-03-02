require 'Proteomics/identifiers'
require 'Proteomics/neighbours'
require 'rbbt/sources/organism'

module Proteomics

  helper :parse_mi do |mi, organism=Organism.default_code("Hsa")|
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
  input :only_pdb, :boolean, "Only consider PDB neighbours", false
  input :just_one, :boolean, "Consider only neighbours from first PDB that contains them", true
  task :mi_neighbours => :tsv do |mis,organism,only_pdb,just_one|

    cpus =  config :cpus, :mi_neighbours, :neighbours, :proteomics, :Proteomics, :default => 2

    annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "PDB", "Neighbours"], :type => :double, :namespace => organism
    annotations.init
    TSV.traverse mis, :cpus => cpus, :bar => self.progress_bar("Mutated Isoform neighbours"), :into => annotations, :type => :array do |mi|

      isoform, residue = parse_mi mi, organism

      next if isoform.nil?

      n = Proteomics.neighbours(isoform, [residue], organism, 5, only_pdb, just_one)

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
end
