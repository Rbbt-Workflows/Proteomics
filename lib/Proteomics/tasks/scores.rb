module Proteomics
  def self.score_for(field, value, all_values = nil)
    value = value.collect{|v| v.split(";")}.flatten
    score = case field
            when "Appris Features"
              if value.include? "firestar"
                2
              else
                1 
              end
            when "UniProt Features"
              relevant = %w(DISULFID DNA_BIND METAL INTRAMEM CROSSLNK MUTAGEN)
              value = value.zip(all_values["UniProt Feature Descriptions"].collect{|v| v.split(";")}.flatten).reject{|v,d| v == "MUTAGEN" and d =~ /no effect/i}.collect{|v,d| v}
              sum = 0
              sum += 1 if (value & relevant).any?
              sum
            when "Sample"
              case 
              when value.length > 10
                3
              when value.length > 5
                2
              when value.length > 1
                1 
              else
                0
              end
            when "UniProt Variant ID"
              case 
              when value.length > 0
                1
              else
                0
              end
            when "Type of Variant"
              if value.include?("Disease")
                2
              elsif value.include?("Unclassified")
                1
              else
                0
              end
            when "Partner Ensembl Protein ID"
              2
            else
              0
            end
    score
  end

  def self.score_mi(values)
    score = 0
    values.zip(values.fields).each do |value, field|
      next if value.empty?
      if field =~ /Neighbour/
        score = score.to_f + (score_for(field.sub('Neighbour ',''), value, values).to_f / 2)
      else
        score = score + score_for(field, value, values)
      end
    end
    score
  end

  dep :wizard
  task :scores => :tsv do 
    wizard = step(:wizard)
    wizard_res = wizard.load

    wizard_res.add_field "Score" do |mi, values|
      Structure.score_mi(values)
    end

    wizard_res.reorder "Mutated Isoform", "Score"
  end
  
  dep :scores
  task :score_summary => :tsv do 
    scores = step(:scores)
    wizard = scores.step(:wizard)
    wizard_res = wizard.load
    scores_res = scores.load


    require 'rbbt/sources/clinvar'
    Workflow.require_workflow "DbNSFP"
    Workflow.require_workflow "Pandrugs"


    clinvar = ClinVar.mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true
    interfaces = wizard.file('interfaces').tsv
    dbNSFP = DbNSFP.job(:annotate, nil, :mutations => wizard_res.keys).run
    cosmic = wizard.file('COSMIC').tsv
    cosmic_neighbours = wizard.file('COSMIC neighbours').tsv
    appris = wizard.file('Appris').tsv
    appris_neighbours = wizard.file('Appris neighbours').tsv
    uniprot = wizard.file('UniProt').tsv
    uniprot_neighbours = wizard.file('UniProt neighbours').tsv
    predictors = %w(SIFT Polyphen2_HDIV Polyphen2_HVAR MutationTaster MutationAssessor FATHMM LRT VEST3 CADD )
    thresholds = %w( <0.05 >0.957,0.453 >0.909,0.447 >0.5 >3.5,1.9 <-1.5 - >0.8 >3.5   )

    fields = %w(Score CV #CS FL MR MUT DT PPI DP)
    tsv = TSV.setup({}, :key_field => "Mutated Isoform", :fields => fields, :type => :list, :namespace => wizard_res.namespace)

    scores_res.each do |mutation, score|
      values = []
      ensp, _sep, change = mutation.partition(":") 
      next unless ensp =~ /^ENSP/

      damage_count = 0
      total_preds = 0
      dvalues = dbNSFP[mutation]
      if dvalues
        predictors.each_with_index do |predictor,i|
          next if predictor == "LRT"
          raw, dscore, converted, rankscore, raw_rankscore, converted_rankscore, p = nil
          threshold = thresholds[i]
          raw = dvalues[predictor + '_raw'] if dvalues.fields.include? predictor + '_raw'
          dscore = dvalues[predictor + '_score'] if dvalues.fields.include? predictor + '_score'
          dscore = nil if String === dscore and dscore.empty?
          dscore = raw if dscore.nil?
          converted = dvalues[predictor + '_converted_score'] if dvalues.fields.include? predictor + '_converted_score'
          rankscore = dvalues[predictor + '_rankscore'] if dvalues.fields.include? predictor + '_rankscore'
          raw_rankscore = dvalues[predictor + '_raw_rankscore'] if dvalues.fields.include? predictor + '_raw_rankscore'
          converted_rankscore = dvalues[predictor + '_converted_rankscore'] if dvalues.fields.include? predictor + '_converted_rankscore'

          if score and threshold != '-'
            p = case threshold
              when /^<(.*)/
                ths = $1.split(",")
                ths.inject(0){|acc,e| acc += 1 if dscore.to_f < e.to_f; acc}.to_f/ths.length
              when /^>(.*)/
                ths = $1.split(",")
                ths.inject(0){|acc,e| acc += 1 if dscore.to_f > e.to_f; acc}.to_f/ths.length
              else
                nil
              end

            damage_count += 1 if p > 0.5
            total_preds +=1

          end
        end
      end

      values << score.flatten.first

      if clinvar.include? mutation and clinvar[mutation] == "Pathogenic"
        values << "Yes"
      else
        values << "No"
      end

      if cosmic.include? mutation
        count = cosmic[mutation]["Sample name"].length
      else
        count = 0
      end

      if cosmic_neighbours.include?(mutation) && ! cosmic_neighbours[mutation]["Sample name"].nil?
        begin
          ncount = cosmic_neighbours[mutation]["Sample name"].collect{|l| l.split(";")}.flatten.uniq.length
        rescue Exception
          raise $!
        end
      else
        ncount = 0
      end

      values << "#{count} (#{ncount})"


      if appris.include? mutation
        count = Misc.zip_fields(appris[mutation]).select{|type,lig| type =~ /firestar/}.length
      else
        count = 0
      end

      if appris_neighbours.include? mutation
        ncount = Misc.zip_fields(appris_neighbours[mutation]).select{|res,type,lig| type =~ /firestar/}.length
      else
        ncount = 0
      end
      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"



      if uniprot.include? mutation
        count = Misc.zip_fields(uniprot[mutation]).select{|feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? mutation
        ncount = Misc.zip_fields(uniprot_neighbours[mutation]).select{|res,feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if uniprot.include? mutation
        count = Misc.zip_fields(uniprot[mutation]).select{|feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? mutation
        ncount = Misc.zip_fields(uniprot_neighbours[mutation]).select{|res,feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if Pandrugs.knowledge_base.subset('gene_drugs', :source => [mutation.protein.gene], :target => :all).filter(:target_marker => 'target').filter(:status => "Approved").length > 0
        values << "Yes"
      else
        values << "No"
      end

      values << interfaces.include?(mutation) ? "Yes" : "No"

      if dbNSFP.include? mutation
        values << "#{damage_count} of 8"
      else
        values << "NA"
      end
    
      tsv[mutation] = values
    end

    tsv
  end
end
