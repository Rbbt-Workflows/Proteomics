require 'rbbt/sources/uniprot'
require 'rbbt/sources/interactome_3d'
require 'Proteomics/identifiers'

module PDB
  def self.isoform_uniprot_pdbs(isoform, organism = Organism.default_organism("Hsa"))
    uniprot = Proteomics.iso2uni(organism)[isoform]
    sequence = Proteomics.iso2seq(organism)[isoform]

    if uniprot.nil?
      []
    else
      UniProt.pdbs(uniprot)
    end
  end
end
