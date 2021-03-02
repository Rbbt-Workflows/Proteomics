require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'rbbt/PDB/isoform'

module Proteomics
  def self.neighbours(isoform, residues, organism = Organism.default_organism("Hsa"), distance = 5, only_pdb = false, just_one = true)
    pdbs = PDB.isoform_uniprot_pdbs(isoform, organism)

    sequence = iso2seq(organism)[isoform]

    tsv = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)

    return tsv if sequence.nil?

    if pdbs
    end

    return tsv if only_pdb

    pdbs.each do |pdb,info|
      begin
        neighbours_in_pdb = neighbours_in_pdb(sequence, residues, pdb, nil, nil, distance)

        next if neighbours_in_pdb.empty?

        neighbours_in_pdb.each do |seq_pos, seq_neigh|
          next if seq_pos.nil?
          tsv[[isoform, seq_pos] * ":"] = [isoform, seq_pos, pdb, seq_neigh * ";"]
        end 

        break if just_one
    
      rescue
        Log.warn "Error processing #{ pdb }: #{$!.message}"
        next
      end
    end

    # TODO Add contiguous also when PDB is present
    residues.each do |p| 
      new = []
      new << p-1 if p > 1
      new << p+1 if p < sequence.length 
      tsv[[isoform,p]*":"] = [isoform, p, nil, new * ";"]
    end unless only_pdb

    tsv
  end
end
