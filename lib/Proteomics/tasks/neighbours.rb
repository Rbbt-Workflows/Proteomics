module Proteomics
  input :distance, :float, "Distance in angstroms", 5
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  def self.neighbour_map(distance, pdb = nil, pdbfile = nil)
    tsv = PDB.close_residues(distance, pdb, pdbfile)

    # Add the reverse contacts
    new = {}
    tsv.each do |p,ns|
      ns.each{|n| new[n] ||= []; new[n] << p }
    end

    new.each{|p,ns| tsv[p] ||= []; tsv[p].concat ns }

    TSV.setup tsv, :key_field => "Residue", :fields => ["Neighbours"], :type => :flat

    tsv
  end
  task :neighbour_map => :tsv 

  input :sequence, :text, "Protein sequence", nil, :required => true
  input :positions, :array, "Positions inside sequence", nil, :required => true
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "Check only a particular chain", nil
  input :distance, :float, "Distance", 5
  def self.neighbours_in_pdb(sequence, positions, pdb = nil, pdbfile = nil, chain = nil, distance = 5)
    sequence.gsub!(/\s/,'')

    neighbours_in_pdb = TSV.setup({}, :key_field => "Sequence position", :fields => ["Neighbours"], :type => :flat)

    positions_in_pdb = Proteomics.sequence_position_in_pdb(sequence, positions, pdb, pdbfile)

    Log.debug "Position in PDB: #{Misc.fingerprint positions_in_pdb}"

    if chain.nil?
      possible_positions = positions_in_pdb.collect.reject{|c,p| p.nil? }
      raise "No matched chain positions" if possible_positions.empty?
      chain =  possible_positions.sort_by{|c,p| p.length}.first.first
    end

    return neighbours_in_pdb if positions_in_pdb.nil? or positions_in_pdb[chain].nil?

    neighbour_map = neighbour_map(distance, pdb, pdbfile)

    alignment_maps = {}
    PDB.pdb_alignment_map(sequence, pdb, pdbfile).each do |chain, map|
      alignment_maps[chain] = map.invert
    end
    alignment_map = alignment_maps[chain]

    positions_in_pdb[chain].each do |position|
      position_in_chain = [chain, position] * ":"
      neigh = neighbour_map[position_in_chain]

      all_neighbours = neigh || []

      position_in_sequence = alignment_map[position.to_i]
      next if position_in_sequence.nil?

      all_neighbours_in_sequence = all_neighbours.collect do |n| 
        c, p = n.split(":")
        alignment_maps[c][p.to_i]
      end.flatten.compact
      neighbours_in_pdb[position_in_sequence] = all_neighbours_in_sequence
    end

    neighbours_in_pdb
  end
  task :neighbours_in_pdb => :tsv

end
