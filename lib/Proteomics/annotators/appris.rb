module Proteomics
  def self.isoform2transcript(organism = Organism.default_code("Hsa"))
    @isoform2transcript ||= Organism.transcripts(organism).index :target => "Ensembl Transcript ID", :fields => "Ensembl Protein ID", :persist => true
  end

  #def self.appris_dataset
  #  @appris_dataset ||= begin
  #                        Rbbt.share.data["appris_data.fir_spa_thu_cra.gen19.v2.tsv"].find(:lib).tsv :key_field => "Ensembl Transcript ID", :persist => true
  #                      end
  #end

  #def self.appris_features(isoform)
  #  transcript = isoform2transcript[isoform]

  #  values = appris_dataset[transcript]
  #  return [] if values.nil?

  #  features = []

  #  values["firestar functional residues"].each do |res|
  #    position, ligand = res.split ":"
  #    features << {:type => "firestar", :start => position.to_i, :end => position.to_i, :description => ligand}
  #  end

  #  values["spade whole domain"].each do |res|
  #    position, pfam_acc = res.split ":"
  #    start, eend = position.split("-")
  #    pfam_acc.sub!(/\.\d+$/, '')
  #    features << {:type => "spade", :start => start.to_i, :end => eend.to_i, :description => pfam_acc}
  #  end

  #  values["thump transmembrane helix"].each do |res|
  #    position, damage = res.split ":"
  #    start, eend = position.split("-")
  #    damage = damage == "1" ? "Damaged" : "Normal"
  #    features << {:type => "thump", :start => start.to_i, :end => eend.to_i, :description => damage}
  #  end

  #  values["crash signal peptide"].each do |res|
  #    start, eend = res.split("-")
  #    features << {:type => "crash", :start => start.to_i, :end => eend.to_i, :description => "Signal peptide"}
  #  end

  #  features
  #end

  #add_annotator("Appris", "Appris Features", "Appris Feature locations", "Appris Feature Descriptions") do |isoform, residue,organism|
  #  features = Proteomics.appris_features(isoform)
  #  next if features.empty?

  #  case residue
  #  when Fixnum
  #    start = eend = residue
  #  when /(\d+):(.*)/
  #    start = $1.to_i
  #    eend = $2.to_i
  #  else
  #    raise "Format of residue not understood: #{residue.inspect}"
  #  end

  #  overlapping = [[],[],[]]
  #  features.select{|info|
  #    (info[:start].to_i <= eend || eend == -1) and info[:end].to_i >= start
  #  }.each{|info|
  #    overlapping[0] << info[:type]
  #    overlapping[1] << [info[:start], info[:end]] * ":"
  #    overlapping[2] << (info[:description] || "").strip.sub(/\.$/,'')
  #  }

  #  next if overlapping.first.empty?
  #  overlapping
  #end

  def self.appris_annotations(organism)
    @@appris_annotations ||= {}
    @@appris_annotations[organism] ||= begin
                                         org, data = organism.split("/")
                                         build = Organism.GRC_build(organism)
                                         protein_annotations = Appris[org][build].protein_annotations.produce.tsv
                                       end
  end

  add_annotator("Appris", "Appris Features", "Appris Feature locations", "Appris Feature Descriptions") do |isoform, residue,organism|
    
    annotations = Proteomics.appris_annotations(organism)
    transcript = Proteomics.isoform2transcript(organism)[isoform]
    next if transcript.nil?
    values = annotations[transcript]
    next if values.nil?
    res = []
    Misc.zip_fields(values).each do |location,feature,value|
      if location.include? "-"
        loc = location.split("-")
        next if loc.first.to_i > residue
        next if loc.last.to_i < residue
      else
        next if location.to_i != residue
      end
      res = [feature, location, value]
    end
    next if res.empty?
    res.extend MultipleResult
    res
  end

end
