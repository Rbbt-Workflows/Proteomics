module Proteomics

  dep :annotate_mi, :database => :placeholder do |jobname,options|
    Proteomics::ANNOTATORS.keys.collect do |database|
      {:inputs => options.merge(:database => database), :jobname => jobname}
    end
  end
  dep :annotate_mi_neighbours, :database => :placeholder do |jobname,options|
    Proteomics::ANNOTATORS.keys.collect do |database|
      {:inputs => options.merge(:database => database), :jobname => jobname}
    end
  end
  dep :mi_interfaces
  task :mi_wizard => :tsv do
    Step.wait_for_jobs dependencies
    dependencies.inject(nil) do |acc,dep|
      tsv = dep.load
      if dep.task_name.to_s.include?("neig")
        tsv.fields = tsv.fields.collect{|f| f == "Neighbour" ? f : "Neighbouring " + f }
      end
      acc = acc.nil? ? tsv : acc.attach(tsv, :complete => true, :zipped => true)
      acc
    end
  end

  dep :annotate_dna, :database => :placeholder, :compute => :produce do |jobname,options|
    Proteomics::ANNOTATORS.keys.collect do |database|
      {:inputs => options.merge(:database => database), :jobname => jobname}
    end
  end
  dep :annotate_dna_neighbours, :database => :placeholder do |jobname,options|
    Proteomics::ANNOTATORS.keys.collect do |database|
      {:inputs => options.merge(:database => database), :jobname => jobname}
    end
  end
  dep :dna_interfaces
  task :dna_wizard => :tsv do
    Step.wait_for_jobs dependencies
    dependencies.inject(nil) do |acc,dep|
      tsv = dep.load
      if dep.task_name.to_s.include?("neig")
        tsv.fields = tsv.fields.collect{|f| f == "Neighbour" ? f : "Neighbouring " + f }
      end
      acc = acc.nil? ? tsv : acc.attach(tsv, :complete => true, :zipped => true)
      acc
    end
  end

  input :mutations, :array, "Mutations (e.g. 18:6237978:G, ENSP00000382976:L257R, L3MBTL4:L257R)"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :wizard => :tsv do |mutations,organism,watson|
    databases = Proteomics::ANNOTATORS.keys 

    raise ParameterException, "No mutations given" if mutations.nil? 
    raise ParameterException, "This wizard is limited to a thousand variants. For larger jobs use standard functions please" if mutations.length > 1000

    log :identifiying
    mutations = mutations.collect{|m| m.sub(/(.*) (.*)/,'\1:\2')}
    type = case mutations.first
           when nil
             raise ParameterException, "No mutations given"
           when /.{1,2}:(\d+):[ACTG+-]+/
             :genomic
           when /ENS.*:\w+\d+\w+/
             :protein
           when /.*:\w+\d+\w+/
             mutations = mutations.collect do |m|
               Proteomics.name2pi(m, organism)
             end.compact
             :protein
           end

    case type
    when :genomic
      all_annotations = Sequence.job(:mutated_isoforms_fast, clean_name, :mutations => mutations, :organism => organism, :principal => true, :watson => watson).run.to_double
      databases.each do |database|
        log database
        annotations = Proteomics.job(:annotate_dna, clean_name, :mutations => mutations, :organism => organism, :database => database, :principal => true, :watson => watson).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end

      log :interfaces
      interfaces = Proteomics.job(:dna_interfaces, clean_name, :mutations => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = Sequence.job(:mutated_isoforms_fast, clean_name, :mutations => mutations, :organism => organism, :principal => true, :watson => watson).run.to_double
      databases.each do |database|
        log database + ' neighbours'
        annotations = Proteomics.job(:annotate_dna_neighbours, clean_name, :mutations => mutations, :organism => organism, :database => database, :principal => true, :watson => watson).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        annotations.rename_field "Residue", database + " residue"
        all_annotations_n = all_annotations_n.attach(annotations)
      end

    when :protein

      all_annotations = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        log database
        annotations = Proteomics.job(:annotate_mi, clean_name, :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end

      log :interfaces
      interfaces = Proteomics.job(:mi_interfaces, clean_name, :mutated_isoforms => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        log database + ' neighbours'
        annotations = Proteomics.job(:annotate_mi_neighbours, clean_name, :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        annotations.rename_field "Residue", database + " residue"
        all_annotations_n = all_annotations_n.attach(annotations)
      end
    end

    all_annotations_n.fields = all_annotations_n.fields.collect{|f| "Neighbour " + f }
    all_annotations.attach(all_annotations_n, :fields => all_annotations_n.fields.reject{|f| f =~ /Mutated Isoform/})

    all_annotations.namespace = organism

    all_annotations
  end

  dep :dna_wizard do |jobname,options|
    if options[:mutations].nil?
      nil
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  dep :mi_wizard do |jobname,options|
    if options[:mutated_isoforms].nil?
      nil
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  task :wizard => :tsv do
    dependencies.first.load
  end
end
