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
        database = dep.recursive_inputs["database"]
        tsv.fields = tsv.fields.collect{|f| f == "Neighbour" ? f + " used for #{database}" : "Neighbouring " + f }
      end
      #  TODO: revise
      #acc = acc.nil? ? tsv : acc.attach(tsv, :complete => true, :zipped => true)
      acc = acc.nil? ? tsv : acc.attach(tsv, :complete => true, :one2one => true)
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
    deps = [step(:mutated_isoforms_fast)] + dependencies
    deps.inject(nil) do |acc,dep|
      tsv = dep.load.to_double
      if dep.task_name.to_s.include?("neig")
        tsv.fields = tsv.fields.collect{|f| f == "Neighbour" ? f + " used for #{database}" : "Neighbouring " + f }
      end
      attach_fields = tsv.fields - ["Genomic Mutation", "Neighbouring Mutated Isoform"]
      acc = acc.nil? ? tsv : acc.attach(tsv, :complete => true, :zipped => true, :fields => attach_fields)
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
