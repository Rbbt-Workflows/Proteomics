require 'rbbt/PDB'

module Proteomics
  input :sequence, :text, "Protein sequence", nil, :required => true
  input :positions, :array, "Positions within protein sequence", nil, :required => true
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  def self.sequence_position_in_pdb(protein_sequence, protein_positions, pdb = nil, pdbfile = nil)
    protein_sequence = protein_sequence.gsub(/\s/,'')
    map = PDB.pdb_alignment_map(protein_sequence, pdb, pdbfile)

    protein_positions = [protein_positions] unless Array === protein_positions
    protein_positions = protein_positions.collect{|p| p.to_i}

    alignments = {}
    map.each do |chain, chain_map|
      alignments[chain] = chain_map.chunked_values_at(protein_positions)
    end

    TSV.setup(alignments, :key_field => "PDB Chain", :fields => protein_positions, :type => :list, :cast => :to_i, :unnamed => true)
  end
  task :sequence_position_in_pdb => :tsv

  input :chain, :string, "PDB chain", nil, :require => true
  input :positions, :array, "Position within PDB chain", nil, :require => true
  input :sequence, :text, "Protein sequence", nil, :require => true
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  def self.pdb_chain_position_in_sequence(chain, positions, protein_sequence, pdb = nil, pdbfile = nil)
    protein_sequence = protein_sequence.gsub(/\s/,'')
    map = PDB.pdb_alignment_map(protein_sequence, pdb, pdbfile)
    chain_map = map[chain]

    map = chain_map.invert

    res = {}
    positions.each do |pos|
      pos = pos.to_i
      res[pos] = map[pos]
    end

    TSV.setup(res, :key_field => "Chain position", :fields => ["Sequence position"], :type => :single, :unnamed => true)
  end
  task :pdb_chain_position_in_sequence => :tsv 


  # In structure/pdb_alignment
  input :sequence, :text, "Protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  def self.pdb_alignment_map(protein_sequence, pdb, pdbfile)
    protein_sequence.gsub!(/\s/,'')


    map = PDB.pdb_alignment_map(protein_sequence, pdb, pdbfile)

    result = TSV.setup({}, :key_field => "Sequence position", :fields => ["Chain:Position in PDB"], :type => :flat)
    map.each do |chain,chain_map|
      chain_map.each do |seq_pos, chain_pos|
        if result[seq_pos].nil?
          result[seq_pos] = [[chain, chain_pos] * ":"]
        else
          result[seq_pos] << [chain, chain_pos] * ":"
        end
      end
    end

    result
  end
  task :pdb_alignment_map => :tsv 
end
