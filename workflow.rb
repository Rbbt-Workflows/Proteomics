require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

Workflow.require_workflow "Appris"

module Proteomics
  extend Workflow


end

require 'Proteomics/tasks/PDB'
require 'Proteomics/tasks/mutated_isoforms'

#require 'rbbt/knowledge_base/Proteomics'
#require 'rbbt/entity/Proteomics'

