require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/sources/clinvar'

Misc.add_libdir if __FILE__ == $0

Workflow.require_workflow "Appris"
Workflow.require_workflow 'COSMIC'
Workflow.require_workflow "InterPro"
Workflow.require_workflow "DbNSFP"
Workflow.require_workflow "Sequence"
#Workflow.require_workflow "Pandrugs"


module Proteomics
  extend Workflow


end

require 'Proteomics/tasks/PDB'
require 'Proteomics/tasks/mutated_isoforms'
require 'Proteomics/tasks/genomic_mutations'
require 'Proteomics/tasks/wizard'
require 'Proteomics/tasks/scores'

Proteomics.tasks.keys.each{|task| Proteomics.export task }
#require 'rbbt/knowledge_base/Proteomics'
#require 'rbbt/entity/Proteomics'

