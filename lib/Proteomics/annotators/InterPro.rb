require 'rbbt/tools/ssw'

module Proteomics

  def self.interpro_protein_domains(organism)
    @interpro_protein_domains ||= InterPro.protein_domains(organism).produce.tsv :persist => true, :unnamed => true
  end

  def self.corrected_interpro_features(uniprot, sequence, organism)
    features = interpro_protein_domains(organism)[uniprot]
    return [] if features.nil? 

    uniprot_sequence = UniProt.sequence(uniprot)

    map = SmithWaterman.alignment_map(uniprot_sequence, sequence)

    corrected_features = []
    Misc.zip_fields(features).each do |code,start,eend|
      corrected_start = map[start.to_i]
      corrected_end = map[eend.to_i]
      next if corrected_start.nil? or corrected_end.nil?
      corrected_info = {}
      corrected_info[:start] = corrected_start
      corrected_info[:end] = corrected_end
      corrected_info[:code] = code
      corrected_features << corrected_info
    end

    corrected_features
  end

  add_annotator("InterPro", "InterPro ID", "Domain range") do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true, :unnamed => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    features =  Misc.insist do
      Persist.persist("Corrected InterPro features", :marshal, :persist => true, :dir => cache(:corrected_interpro_features), :other => {:uniprot => uniprot, :sequence => sequence, :organism => organism}) do
        Proteomics.corrected_interpro_features(uniprot, sequence, organism)
      end
    end
    next if features.empty?

    case residue
    when Integer
      start = eend = residue
    when /(\d+):(.*)/
      start = $1.to_i
      eend = $2.to_i
    else
      raise "Format of residue not understood: #{residue.inspect}"
    end


    overlapping = [[],[]]
    features.select{|info|
      (info[:start].to_i <= eend || eend == -1) and info[:end].to_i >= start
    }.each{|info|
      overlapping[0] << info[:code]
      overlapping[1] << [info[:start], info[:end]] * ":"
    }

    next if overlapping.first.empty?
    overlapping
  end
end
