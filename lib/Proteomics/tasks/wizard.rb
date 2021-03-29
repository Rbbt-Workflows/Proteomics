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
