require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'rbbt/PDB'
require 'rbbt/PDB/isoform'

module Proteomics

  def self.interface_neighbours_i3d(isoform, positions, organism = Organism.default_code("Hsa"), distance = 8)
    tsv = TSV.setup({}, :key_field => "Isoform", :fields => ["Position", "Partner (Ensembl Protein ID)", "PDB", "Partner residues"], :type => :double)

    pdbs = PDB.isoform_i3d_interaction_pdb(isoform, organism)
    sequence = Proteomics.iso2seq(organism)[isoform]

    TSV.traverse pdbs do |filename,values|
      Misc.zip_fields(values).each do |type, _isoform, partner, chain, partner_chain, url|

        partner_positions = []
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

          partner_positions << partner_position
        end

        next if partner_positions.empty?

        map = PDB.pdb_alignment_map(sequence, url)

        partner_map = PDB.pdb_alignment_map(Proteomics.iso2seq(organism)[partner], url)

        residue = map.invert[position]
        partner_residues = partner_map.invert.chunked_values_at(partner_positions)

        tsv.zip_new isoform, [residue, partner, filename, partner_residues]
      end
    end
    tsv
  end


  def self.neighbours(isoform, residues, organism = Organism.default_organism("Hsa"), distance = 5, only_pdb = false, just_one = true)
    pdbs = PDB.isoform_uniprot_pdbs(isoform, organism)

    sequence = iso2seq(organism)[isoform]

    tsv = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)

    return tsv if sequence.nil?

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
      tsv[[isoform, p]*":"] ||= [isoform, p, nil, new * ";"]
    end unless only_pdb

    tsv
  end

  I3D_PROTEINS = Interactome3d.proteins_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS = Interactome3d.interactions_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS_REVERSE = Interactome3d.interactions_tsv.tsv :merge => true, :key_field => "PROT2", :zipped => true, :unnamed => true, :persist => true

  def self.___interface_neighbours_i3d(protein, positions, organism = Organism.default_code("Hsa"), distance = 8)
    tsv = TSV.setup({}, :key_field => "Isoform", :fields => ["Position", "Partner Ensembl Protein ID", "PDB", "Partner residues"], :type => :double)







    uniprot = iso2uni(organism)[protein]
    return tsv if uniprot.nil?
    sequence = iso2seq(organism)[protein]
    return tsv if sequence.nil?

    uni2iso = uni2iso(organism)

    forward_positions = ["PDB_ID", "PROT2", "CHAIN1", "CHAIN2", "FILENAME"].collect{|f| I3D_INTERACTIONS.identify_field f}
    reverse_positions = ["PDB_ID", "PROT1", "CHAIN2", "CHAIN1", "FILENAME"].collect{|f| I3D_INTERACTIONS_REVERSE.identify_field f}

    {:forward => I3D_INTERACTIONS,
      :reverse => I3D_INTERACTIONS_REVERSE}.each do |db_direction, db|

      next unless db.include? uniprot
      values_list = db[uniprot]
      seen_pbs = []
      Misc.zip_fields(values_list).each do |values|
        if db_direction == :forward
          pdb, partner, orig_chain, orig_partner_chain, filename = values.values_at *forward_positions
          chain, partner_chain = "A", "B"
        else
          pdb, partner, orig_chain, orig_partner_chain, filename = values.values_at *reverse_positions
          chain, partner_chain = "B", "A"
        end

        next if chain.strip.empty? or partner_chain.strip.empty?

        if uniprot == partner
          partner_ensembl = protein
        else
          partner_ensembl =  uni2iso[partner]
        end

        if partner_ensembl.nil?
          Log.warn "Could not translate partner to Ensembl: #{ partner }"
          next
        end

        partner_sequence = iso2seq[partner_ensembl]
        if partner_sequence.nil?
          Log.warn "Could get partner sequence: #{ partner } (#{partner_ensembl})"
          next
        end


        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=#{ type }&pdb=#{ filename }"
        Log.debug "Processing: #{ url }"

        positions_in_pdb = sequence_position_in_pdb(sequence, positions, url, nil)[chain]
        next if positions_in_pdb.nil? or positions_in_pdb.empty?

        map = neighbour_map_job url, nil, distance
        #map.unnamed = true

        next if map.nil? or map.empty?

        positions_in_pdb.zip(positions).each do |pdb_position, position|
          code = [chain,pdb_position]*":"
          begin
            neighbours = map[code]
          rescue
            Log.exception $!
            next
          end
          next if neighbours.nil? or neighbours.empty?
          partner_neighbours = neighbours.select{|c| c.split(":").first == partner_chain }.collect{|c| c.split(":").last}
          next if partner_neighbours.empty?
          partner_neighbours_in_sequence = pdb_chain_position_in_sequence(url, nil, partner_chain, partner_neighbours, partner_sequence).values.compact.flatten
          next if partner_neighbours_in_sequence.empty?
          tsv.zip_new(protein, [position, partner_ensembl, url, partner_neighbours_in_sequence * ";"])
        end

        if not seen_pbs.include? pdb
          next if orig_chain == orig_partner_chain
          seen_pbs << pdb
          url = pdb
          Log.debug "Processing: #{ url }"

          positions_in_pdb = sequence_position_in_pdb(sequence, positions, url, nil)[orig_chain]
          next if positions_in_pdb.nil? or positions_in_pdb.empty?

          map = neighbour_map_job url, nil, distance
          #map.unnamed = true

          next if map.nil? or map.empty?

          positions_in_pdb.zip(positions).each do |pdb_position, position|
            code = [orig_chain,pdb_position]*":"
            begin
              neighbours = map[code]
            rescue
              Log.exception $!
              next
            end
            next if neighbours.nil? or neighbours.empty?
            partner_neighbours = neighbours.select{|c| c.split(":").first == orig_partner_chain }.collect{|c| c.split(":").last}
            next if partner_neighbours.empty?
            partner_neighbours_in_sequence = pdb_chain_position_in_sequence(url, nil, orig_partner_chain, partner_neighbours, partner_sequence).values.compact.flatten
            next if partner_neighbours_in_sequence.empty?
            tsv.zip_new(protein, [position, partner_ensembl, url, partner_neighbours_in_sequence * ";"])
          end
        end
      end
    end

    tsv
  end
end
