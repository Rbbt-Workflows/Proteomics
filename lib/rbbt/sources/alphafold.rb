require 'rbbt-util'

module AlphaFold
  BASE_URL="https://alphafold.ebi.ac.uk/api"
  def self.json(uniprot)
    begin
      JSON.parse(Open.read("#{BASE_URL}/prediction/#{uniprot}"))
    rescue
      Log.warn "Could not find entry in AlphaFold: #{uniprot}"
      []
    end
  end

  def self.pdbs(uniprot)
    info = json(uniprot)
    return [] if info.nil? || info.empty?
    
    info.collect{|i| i['pdbUrl'] }
  end
end
