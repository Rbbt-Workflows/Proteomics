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
      [Proteomics.score_mi(values).to_i]
    end

    if wizard_res.key_field == "Genomic Mutation"
      wizard_res.slice ["Mutated Isoform", "Score"]
    else
      wizard_res.slice ["Score"]
    end
  end
  
  dep :scores
  task :score_summary => :tsv do 
    scores = step(:scores)
    wizard = scores.step(:wizard)

    organism = recursive_inputs[:organism]

    wizard_res = wizard.load
    scores_res = scores.load

    by_dna = wizard_res.key_field == "Genomic Mutation"

    #if by_dna
    #  wizard_res = wizard_res.reorder "Mutated Isoform", wizard_res.fields
    #  scores_res = scores_res.reorder "Mutated Isoform", scores_res.fields
    #end

    if by_dna
      mis = wizard_res.column("Mutated Isoform").values.flatten.uniq
    else
      mis = wizard_res.keys
    end

    build = Organism.hg_build(organism)

    clinvar = ClinVar[build].mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true, :monitor => true

    if by_dna
      interfaces = wizard.step('dna_interfaces').load
      cosmic = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna" && dep.inputs[:database] == 'COSMIC' }.first.load
      cosmic_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna_neighbours" && dep.inputs[:database] == 'COSMIC' }.first.load
      appris = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna" && dep.inputs[:database] == 'Appris' }.first.load
      appris_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna_neighbours" && dep.inputs[:database] == 'Appris' }.first.load
      uniprot = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna" && dep.inputs[:database] == 'UniProt' }.first.load
      uniprot_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_dna_neighbours" && dep.inputs[:database] == 'UniProt' }.first.load
    else
      interfaces = wizard.step('mi_interfaces').load
      cosmic = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi" && dep.inputs[:database] == 'COSMIC' }.first.load
      cosmic_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi_neighbours" && dep.inputs[:database] == 'COSMIC' }.first.load
      appris = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi" && dep.inputs[:database] == 'Appris' }.first.load
      appris_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi_neighbours" && dep.inputs[:database] == 'Appris' }.first.load
      uniprot = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi" && dep.inputs[:database] == 'UniProt' }.first.load
      uniprot_neighbours = wizard.rec_dependencies.select{|dep| dep.task_name.to_s == "annotate_mi_neighbours" && dep.inputs[:database] == 'UniProt' }.first.load
    end

    mis_ensp = mis.collect{|mi|
      mi =~ /ENSP/ ? mi : Proteomics.name2pi(mi)
    }

    dbNSFP = DbNSFP.job(:annotate, nil, :mutations => mis_ensp).run
    dbNSFP_pred = DbNSFP.job(:predict, nil, :mutations => mis_ensp).run

    Open.write(file("dbNSFP_pred.tsv"), dbNSFP_pred.to_s)

    predictors = %w(SIFT Polyphen2_HDIV Polyphen2_HVAR MutationTaster MutationAssessor FATHMM LRT VEST3 CADD )
    thresholds = %w( <0.05 >0.957,0.453 >0.909,0.447 >0.5 >3.5,1.9 <-1.5 - >0.8 >3.5   )

    fields = ["Score", "ClinVar", "COSMIC", "Appris", "UniProt MOD_RES", "UniProt MUTAGEN", "PanDrugs", "Affected PPI", "Damage Predictions"]
    tsv = TSV.setup({}, :key_field => "Mutated Isoform", :fields => fields, :type => :list, :namespace => wizard_res.namespace)

    scores_res.each do |key, values|
      if by_dna
        mi, score = values
      else
        score = values.first
        mi = key
      end
      
      mi = mi.first if Array ===  mi

      values = []
      ensp_mi = mi =~ /ENSP/ ? mi : Proteomics.name2pi(mi)
      next unless ensp_mi =~ /^ENSP/


      damage_count = 0
      total_preds = 0
      dvalues = dbNSFP_pred[ensp_mi]

      if dvalues
        damage_count = dvalues.select{|v| %(D H).include? v.to_s}.length
        total_preds = dvalues.select{|v| v.to_s != ""}.length
        #predictors.each_with_index do |predictor,i|
        #  next if predictor == "LRT"
        #  raw, dscore, converted, rankscore, raw_rankscore, converted_rankscore, p = nil
        #  threshold = thresholds[i]
        #  raw = dvalues[predictor + '_raw'] if dvalues.fields.include? predictor + '_raw'
        #  dscore = dvalues[predictor + '_score'] if dvalues.fields.include? predictor + '_score'
        #  dscore = nil if String === dscore and dscore.empty?
        #  dscore = raw if dscore.nil?
        #  converted = dvalues[predictor + '_converted_score'] if dvalues.fields.include? predictor + '_converted_score'
        #  rankscore = dvalues[predictor + '_rankscore'] if dvalues.fields.include? predictor + '_rankscore'
        #  raw_rankscore = dvalues[predictor + '_raw_rankscore'] if dvalues.fields.include? predictor + '_raw_rankscore'
        #  converted_rankscore = dvalues[predictor + '_converted_rankscore'] if dvalues.fields.include? predictor + '_converted_rankscore'

        #  if score and threshold != '-'
        #    p = case threshold
        #      when /^<(.*)/
        #        ths = $1.split(",")
        #        ths.inject(0){|acc,e| acc += 1 if dscore.to_f < e.to_f; acc}.to_f/ths.length
        #      when /^>(.*)/
        #        ths = $1.split(",")
        #        ths.inject(0){|acc,e| acc += 1 if dscore.to_f > e.to_f; acc}.to_f/ths.length
        #      else
        #        nil
        #      end

        #    damage_count += 1 if p > 0.5
        #    total_preds +=1
        #  end
        #end
      end

      values << score.flatten.first

      if clinvar.include? mi and clinvar[mi] == "Pathogenic"
        values << "Yes"
      else
        values << "No"
      end

      if cosmic.include? key
        count = cosmic[key]["Sample name"].length
      else
        count = 0
      end

      if cosmic_neighbours.include?(key) && ! cosmic_neighbours[key]["Sample name"].nil?
        begin
          ncount = cosmic_neighbours[key]["Sample name"].collect{|l| l.split(";")}.flatten.uniq.length
        rescue Exception
          raise $!
        end
      else
        ncount = 0
      end

      values << "#{count} (#{ncount})"

      if appris.include? key
        count = Misc.zip_fields(appris[key]).select{|type,lig| type !~ /wrong|damaged/}.length
      else
        count = 0
      end

      if appris_neighbours.include? key
        ncount = Misc.zip_fields(appris_neighbours[key]).select{|res,type,range,lig| type !~ /wrong|damaged/}.collect{|res,type,range,lig| [type,range,lig]}.uniq.length
      else
        ncount = 0
      end
      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if uniprot.include? key
        count = Misc.zip_fields(uniprot[key]).select{|feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? key
        ncount = Misc.zip_fields(uniprot_neighbours[key]).select{|res,feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if uniprot.include? key
        count = Misc.zip_fields(uniprot[key]).select{|feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? key
        ncount = Misc.zip_fields(uniprot_neighbours[key]).select{|res,feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      begin
        Workflow.require_workflow "Pandrugs"
        ensp = mi.partition(":").first
        gene = Proteomics.ensp2ensg[ensp]
        if Pandrugs.knowledge_base.subset('gene_drugs', :source => [gene], :target => :all).filter(:target_marker => 'target').filter(:status => "Approved").length > 0
          values << "Yes"
        else
          values << "No"
        end
      rescue
        Log.exception $!
        Log.low "No Pandrugs"
        values << "No"
      end

      values << interfaces.include?(mi) ? "Yes" : "No"

      if dbNSFP.include? ensp_mi
        values << "#{damage_count} of #{total_preds}"
      else
        values << "NA"
      end
    
      tsv[mi] = values
    end

    tsv
  end
end
