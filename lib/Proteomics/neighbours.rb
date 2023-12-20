require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'rbbt/sources/alphafold'
require 'rbbt/PDB'
require 'rbbt/PDB/isoform'

module Proteomics

  def self.interface_neighbours_i3d(isoform, residues, organism = Organism.default_code("Hsa"), distance = 8)
    tsv = TSV.setup({}, :key_field => "Isoform", :fields => ["Position", "Partner (Ensembl Protein ID)", "PDB", "Partner residues"], :type => :double)

    pdbs = PDB.isoform_i3d_interaction_pdb(isoform, organism)
    sequence = Proteomics.iso2seq(organism)[isoform]

    return tsv if sequence.nil?

    partners_info = {}
    TSV.traverse pdbs do |filename,values|
      Misc.zip_fields(values).each do |_isoform, partner, chain, partner_chain, type, url|

        partner_positions = {}
        PDB.interface_residues(url, nil, distance).each do |r1, r2|
          c1, p1 = r1.split(":")
          c2, p2 = r2.split(":")

          position = case chain
                     when c1
                       p1
                     when c2
                       p2
                     else
                       next
                     end

          partner_position = case partner_chain
                             when c1
                               p1
                             when c2
                               p2
                             else
                               next
                             end

          partner_positions[[chain, position]] ||= []
          partner_positions[[chain, position]] << [partner_chain, partner_position]
        end


        next if partner_positions.empty?


        map = PDB.pdb_alignment_map(sequence, url)

        parner_sequence = Proteomics.iso2seq(organism)[partner]
        partner_map = PDB.pdb_alignment_map(parner_sequence, url)

        partner_positions.each do |cposition, partner_cposition_list|
          chain, position = cposition
          next if map[chain].nil?
          residue = map[chain].invert[position.to_i]
          next unless residues.include? residue
          partner_residues = partner_cposition_list.collect do |partner_chain, partner_position|
            next if partner_map[partner_chain].nil?
            partner_map[partner_chain].invert[partner_position.to_i]
          end.flatten.compact
          next if partner_residues.empty?
          partners_info[partner] ||= {}
          partners_info[partner][residue] ||= []
          partners_info[partner][residue] << [filename, partner_residues]
        end
      end
    end
    partners_info.each do |partner,resinfo|
      resinfo.each do |residue,list|
        files, partner_residues = Misc.zip_fields(list)
        tsv.zip_new isoform, [residue, partner, files.uniq * ";", partner_residues.uniq * ";"]
      end
    end
    tsv
  end


  def self.neighbours(isoform, residues, organism = Organism.default_organism("Hsa"), distance = 5, only_pdb = false, just_one = true, use_i3d = true, use_alphafold = true)
    tsv = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :double)

    pdbs = PDB.isoform_uniprot_pdbs(isoform, organism).keys

    if use_alphafold
      uniprot = Proteomics.iso2uni(organism)[isoform]
      pdbs += pdbs + AlphaFold.pdbs(uniprot) if uniprot
    end

    pdbs += pdbs + PDB.isoform_i3d_pdbs(isoform, organism).column("URL").values.flatten if use_i3d
    
    sequence = iso2seq(organism)[isoform]

    return tsv if sequence.nil?

    pdbs.each do |pdb,info|
      begin
        neighbours_in_pdb = neighbours_in_pdb(sequence, residues, pdb, nil, nil, distance)

        next if neighbours_in_pdb.empty?

        neighbours_in_pdb.each do |seq_pos, seq_neigh|
          next if seq_pos.nil?
          tsv.zip_new([isoform, seq_pos] * ":", [isoform, seq_pos, pdb, seq_neigh * ";"])
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
      tsv.zip_new([isoform, p]*":", [isoform, p, nil, new * ";"])
    end unless only_pdb

    tsv
  end

end
