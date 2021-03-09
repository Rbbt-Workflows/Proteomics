require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/interactome_3d'
require 'Proteomics/identifiers'

module PDB
  def self.isoform_uniprot_pdbs(isoform, organism = Organism.default_code("Hsa"))
    uniprot = Proteomics.iso2uni(organism)[isoform]
    return {} if uniprot.nil?

    UniProt.pdbs(uniprot)
  end

  def self.isoform_i3d_pdbs(isoform, organism = Organism.default_code("Hsa"))
    uniprot = Proteomics.iso2uni(organism)[isoform]

    pdbs = TSV.setup({}, "PDB~Chain,Type,URL#:type=:double")

    return pdbs if uniprot.nil?

    proteins = Interactome3d.proteins_tsv.tsv :merge => true, :persist => true

    values = proteins[uniprot]

    return pdbs if values.nil?

    Misc.zip_fields(values).each do |list|
      info = Misc.zip2hash(values.fields, list)
      chain, filename, type = info.values_at "CHAIN", "FILENAME", "TYPE"

      type = filename =~ /EXP/ ? :pdb : :model
      url = "https://interactome3d.irbbarcelona.org//pdb.php?dataset=human&type1=proteins&type2=#{type}&pdb=#{filename}"

      pdbs.zip_new(filename, [chain, type, url])
    end

    pdbs
  end

  def self.isoform_i3d_interaction_pdb(isoform, organism = Organism.default_code("Hsa"))
    uniprot = Proteomics.iso2uni(organism)[isoform]

    pdbs = TSV.setup({}, "PDB~Protein (Ensembl Protein ID),Partner (Ensembl Protein ID),Chain,Partner Chain,Type,URL#:type=:double")

    return pdbs if uniprot.nil?

    interactions = @@interactions ||= Interactome3d.interactions_tsv.tsv :merge => true, :persist => true
    interactions_reverse = @@interactions_reverse ||= Interactome3d.interactions_tsv.tsv :merge => true, :key_field => "PROT2", :zipped => true, :persist => true

    if values = interactions[uniprot]
      prot1 = uniprot
      Misc.zip_fields(values).each do |list|
        info = Misc.zip2hash(values.fields, list)
        prot2, chain1, chain2, filename, type = info.values_at "PROT2", "CHAIN1", "CHAIN2", "FILENAME", "TYPE"
        
        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=#{ type }&pdb=#{ filename }"

        if prot2 == prot1
          partner = isoform
        else
          partner = Proteomics.uni2iso(organism)[prot2]
        end
        pdbs.zip_new(filename, [isoform, partner, chain1, chain2, type, url])
      end
    end

    if values = interactions_reverse[uniprot]
      prot2 = uniprot
      Misc.zip_fields(values).each do |list|
        info = Misc.zip2hash(values.fields, list)
        prot1, chain1, chain2, filename, type = info.values_at "PROT1", "CHAIN1", "CHAIN2", "FILENAME", "TYPE"
        
        if prot2 == prot1
          next
          partner = isoform
        else
          partner = Proteomics.uni2iso(organism)[prot1]
        end

        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=#{ type }&pdb=#{ filename }"

        pdbs.zip_new(filename, [isoform, partner, chain2, chain1, type, url])
      end
    end

    pdbs
  end
end
