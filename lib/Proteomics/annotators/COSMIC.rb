Workflow.require_workflow 'COSMIC'

module Proteomics

  def self.COSMIC_residues
    @@COSMIC_residues ||= Persist.persist_tsv(nil, "COSMIC::residues", {}, :persist => true, :serializer => :list, :dir => Annotator.cache(:cosmic_residues)) do |data|
                           isoform_residue_mutations = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Genomic Mutations"], :type => :flat)

                           db = COSMIC.knowledge_base.get_database('mutation_isoforms')

                           db.monitor = {:desc => "Processing COSMIC", :step => 10000}
                           db.with_unnamed do
                            db.through  do |mutation, mis|
                              mis.flatten.each do |mi|
                                next unless mi =~ /(ENSP\d+):([A-Z])(\d+)([A-Z]|Indel|Frameshift)$/ and $2 != $4
                                residue = $3.to_i
                                protein = $1
                                isoform_residue_mutations[[protein, residue] * ":"] ||= []
                                isoform_residue_mutations[[protein, residue] * ":"] << mutation
                              end
                            end
                           end

                           data.merge! isoform_residue_mutations
                           isoform_residue_mutations.annotate data

                           data
                         end
  end

  def self.COSMIC_mutation_annotations
    @@COSMIC_mutation_annotations ||= begin
                                       fields = [
                                         'Sample name',
                                         'Primary site',
                                         'Site subtype 1',
                                         'Site subtype 2',
                                         'Site subtype 3',
                                         'Primary histology',
                                         'Histology subtype 1',
                                         'Histology subtype 2',
                                         'Histology subtype 3',
                                         'Pubmed_PMID',
                                       ]
                                       Persist.persist_tsv(COSMIC.mutations, nil, { :key_field => "Genomic Mutations", :fields => fields}, :persist => true ) do |data|
                                         #COSMIC.mutations.tsv :key_field => "Genomic Mutation", :fields => fields, :persist => true, :unamed => true, :type => :double, :zipped => true, :monitor => true
                                         reorder = TSV.reorder_stream_tsv COSMIC.mutations, "Genomic Mutation", fields, true, "Reorder COSMIC mutations"
                                         collapsed = TSV.collapse_stream reorder
                                         TSV.open(collapsed.stream, :key_field => "Genomic Mutation", :persist => true, :persist_data => data, :unamed => true, :type => :double, :monitor => true)
                                       end
                                     end
  end

  def self.COSMIC_resistance_mutations
    @@COSMIC_resistance_mutations ||= begin
                                        fix_change = lambda{|line|
                                          if line[0] == "#"
                                            line
                                          else
                                            mi, *rest = line.chomp.split("\t", -1)
                                            re = mi.match(/^(.*):([A-Z*?])(\d+)([A-Z*?]+)$/)
                                            raise TSV::Parser::SKIP_LINE if re.nil?

                                            isoform = re[1]
                                            residue = re[3]
                                            key = [isoform, residue] * ":"

                                            rest.unshift key
                                            rest * "\t"
                                          end
                                        }
                                        COSMIC.mi_drug_resistance.tsv :fix => fix_change, :merge => true, :persist => true, :unnamed => true
                                      end
  end


  def self.COSMIC_complete_cna
    @@COSMIC_complete_cna ||= begin
                                Persist.persist_tsv(COSMIC.completeCNA, nil, {:key_field => "Sample name:Ensembl Protein ID", :fields => ["CNA"]}, {:persist => true, :serializer => :single}) do |data|
                                  organism = "Hsa/feb2014"
                                  enst2ensp = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :merge => true, :type => :flat, :persist => true
                                  TSV.traverse COSMIC.completeCNA, :bar => true, :type => :double, :into => data do |sample, values|
                                    transcripts, cnas = values
                                    proteins = enst2ensp.chunked_values_at transcripts
                                    res = []
                                    proteins.zip(cnas).each do |ensp,cna|
                                      key = [sample, ensp] * ":"
                                      res << [key, cna]
                                    end
                                    res.extend MultipleResult
                                    res
                                  end
                                end
                              end
  end

  def self.COSMIC_complete_gene_expression
    @@COSMIC_complete_gene_expression ||= begin
                                Persist.persist_tsv(COSMIC.geneExpression, nil, {:key_field => "Sample name:Ensembl Protein ID", :fields => ["Regulation"]}, {:persist => true, :serializer => :single}) do |data|
                                  organism = "Hsa/feb2014"
                                  enst2ensp = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :merge => true, :type => :flat, :persist => true

                                  TSV.traverse COSMIC.geneExpression, :bar => true, :type => :double, :cpus => 10 do |sample, values|
                                    transcripts, exprs = values
                                    proteins = enst2ensp.chunked_values_at transcripts
                                    res = []
                                    proteins.zip(exprs).each do |ensp,expr|
                                      key = [sample, ensp] * ":"
                                      data[key] = expr
                                    end
                                  end
                                end
                              end
  end

  add_annotator("COSMIC", "Genomic Mutation", 'Sample name', 'Primary site', 'Site subtype 1', 'Site subtype 2', 'Site subtype 3', 'Primary histology', 'Histology subtype 1', 'Histology subtype 2', 'Histology subtype 3', "PMID") do |isoform, residue, organism|

    @cosmic_residue_mutations ||= Proteomics.COSMIC_residues
    @cosmic_mutation_annotations ||= Proteomics.COSMIC_mutation_annotations

    isoform_residue = isoform + ":" << residue.to_s
    mutations = @cosmic_residue_mutations[isoform_residue]

    next if mutations.nil?
    mutations.uniq!
    next if mutations.empty?
    tmp = {}
    mutations.each do |mutation|
      annot = @cosmic_mutation_annotations[mutation]
      next if annot.nil?
      Misc.zip_fields(annot).each do |a|
        sample, *rest = a
        sample_isoform = sample + ":" + isoform
        next if tmp.include? sample
        tmp[sample] = rest
      end
    end
    annot = tmp.collect{|p| p.flatten}
    annot = Misc.zip_fields(annot)
    annot.unshift mutations
    annot
  end

  add_annotator("COSMIC_resistance_mutations", "Drug", "Sample", "PMID", "Zygosity") do |isoform,residue,organism|
    db = Proteomics.COSMIC_resistance_mutations
    isoform_residue = [isoform, residue] * ":"
    db[isoform_residue]
  end
end
