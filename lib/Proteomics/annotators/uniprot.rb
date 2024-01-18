require 'rbbt/sources/uniprot'
require 'rbbt/tools/ssw'

module Proteomics
  def self.uniprot_sequence_map(uniprot, sequence)
    uniprot_sequence = UniProt.sequence(uniprot)
    SmithWaterman.alignment_map(uniprot_sequence, sequence)
  end

  def self.corrected_uniprot_features(uniprot, sequence)
    Persist.persist("Corrected UniProt features", :marshal,  :dir => Proteomics::Annotator.cache(:corrected_uniprot_features), :other => {:uniprot => uniprot, :sequence => sequence}) do
      uniprot_sequence = UniProt.sequence(uniprot)

      map = uniprot_sequence_map(uniprot, sequence)

      features = UniProt.features(uniprot)
      corrected_features = []
      features.each do |info|
        corrected_start = map[info[:start]]
        corrected_end = map[info[:end]]
        next if corrected_start.nil? or corrected_end.nil?
        corrected_info = info.dup
        corrected_info[:start] = corrected_start
        corrected_info[:end] = corrected_end
        corrected_features << corrected_info
      end

      corrected_features
    end
  end


  def self.UniProt_residues
    @UniProt_residues ||= Persist.persist_tsv(UniProt.annotated_variants, "UniProt::residues", {}, :persist => true, :serializer => :list, :dir => Rbbt.var.persistence.find(:lib)) do |data|
                           isoform_residue_mutations = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["UniProt Variant ID"], :type => :flat)

                           uni2ensp = Organism.protein_identifiers(Organism.default_code("Hsa")).tsv :fields => ["Ensembl Protein ID"], :key_field => "UniProt/SwissProt Accession", :persist => true, :type => :flat, :merge => true, :unnamed => true
                           ensp2sequence = Organism.protein_sequence(Organism.default_code("Hsa")).tsv :persist => true, :unnamed => true

                           db = UniProt.annotated_variants.tsv(:fields => ["Amino Acid Mutation", "UniProt Variant ID"], :persist => true, :type => :double, :unnamed => true)
                           db.monitor = {:desc => "Processing UniProt", :step => 1000}

                           db.with_unnamed do
                             db.through do |uniprot,values|
                               Log.debug Log.color :red, uniprot
                               begin
                                 ensps = uni2ensp[uniprot]
                                 raise "No translation to Ensembl: #{ uniprot }" if ensps.nil? or ensps.empty?
                                 found = false
                                 ensps.each do |ensp|
                                   Log.debug Log.color :blue, ensp
                                   begin
                                     ensp_sequence = ensp2sequence[ensp]
                                     raise "No sequence: #{ ensp } " if ensp_sequence.nil?
                                     uniprot_sequence = UniProt.sequence(uniprot)
                                     map = SmithWaterman.alignment_map(uniprot_sequence, ensp_sequence)

                                     Misc.zip_fields(values).each do |change,vid|
                                       match = change.match(/^([A-Z])(\d+)([A-Z])$/)
                                       raise "Unknown change: #{ ensp } #{change}" if match.nil?
                                       ref, _pos, mut = match.values_at 1,2,3
                                       pos = map[_pos.to_i]
                                       raise "Unmapped position: #{ ensp } #{_pos}" if pos.nil?
                                       isoform_residue_mutations[[ensp,pos]*":"] ||= []
                                       isoform_residue_mutations[[ensp,pos]*":"] << vid
                                       found = true
                                     end
                                   rescue
                                     Log.debug $!.message
                                     next
                                   end
                                 end
                                 raise "No suitable mapings for #{ uniprot }" unless found
                               rescue
                                 Log.warn $!.message
                                 next
                               end
                             end
                           end

                           Log.info "Merging data"
                           data.merge! isoform_residue_mutations
                           Log.info "Annotating database"
                           isoform_residue_mutations.annotate data

                           data
    end
  end

  def self.UniProt_mutation_annotations
    @UniProt_mutation_annotations ||= begin
                                        fields = [
                                          'Amino Acid Mutation',
                                          'Type of Variant',
                                          'SNP ID',
                                          'Disease'
                                        ]
                                       #UniProt.annotated_variants.tsv(:key_field => "UniProt Variant ID", :fields => fields, :persist => true, :type => :double, :unnamed => true, :one2one => true, :zipped => true).to_list
                                       UniProt.annotated_variants.tsv(:key_field => "UniProt Variant ID", :fields => fields, :persist => true, :type => :double, :unnamed => true, :one2one => true).to_list
                                       raise
                                     end
  end
  
  Proteomics.add_annotator("UniProt", "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions") do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true, :unnamed => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    _other = {:uniprot => uniprot, :sequence => sequence}
    features = Proteomics.corrected_uniprot_features(uniprot, sequence)

    next if features.empty?

    overlapping = [[],[],[]]

    case residue
    when Integer
      start = eend = residue
    when /(\d+):(.*)/
      start = $1.to_i
      eend = $2.to_i
    else
      raise "Format of residue not understood: #{residue.inspect}"
    end

    features.select{|info|
      case info[:type]
      when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
        false
      when "DISULFID", "CROSSLNK", "VARIANT"
        ([info[:start], info[:end]] & [start, eend]).any?
      else
        (info[:start].to_i <= eend || eend == -1) and info[:end].to_i >= start
      end
    }.each{|info|
      description = (info[:description] || "").strip.sub(/\.$/,'')
      overlapping[0] << info[:type]
      overlapping[1] << [info[:start], info[:end]] * ":"
      overlapping[2] << description.gsub('|','-').gsub(';','-')
    }

    next if overlapping.first.empty?

    overlapping
  end

  add_annotator("variants", "UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease") do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    @annotations ||= Proteomics.UniProt_mutation_annotations
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true, :unnamed => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    features =  Proteomics.corrected_uniprot_features(uniprot, sequence)

    next if features.empty?

    overlapping = [[],[],[],[]]

    features.select{|info|
      info[:type] == "VARIANT" and info[:start] == residue
    }.each{|info|
      if m = info[:description].match(/(VAR_\d+)/)
        id = m[1]
        next unless @annotations.include?(id)
        overlapping[0] << id
        annots = @annotations[id]
        overlapping[1] << annots[2]
        overlapping[2] << annots[1]
        overlapping[3] << annots[3]
      end
    }

    next if overlapping.first.empty?

    overlapping
  end

end
