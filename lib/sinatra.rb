require 'rbbt/rest/web_tool'
include Sinatra::RbbtToolHelper

require 'rbbt-util'
require 'rbbt/workflow'

#require 'rbbt/entity/structure'

Workflow.require_workflow "Genomics"
require 'rbbt/entity/gene'
require 'rbbt/entity/protein'
require 'rbbt/entity/mutated_isoform'
require 'rbbt/entity/InterPro'
require 'rbbt/sources/pfam'

[Gene, Protein, MutatedIsoform, InterProDomain, PfamDomain].each do |mod|
  mod.instance_eval do
    extend Entity
    include Entity::REST
  end
end

MutatedIsoform.instance_eval do
  property :normalized => :single do
    gene, change = self.split(":")
    protein = Proteomics.gene2isoform(gene, organism)
    MutatedIsoform.setup([protein, change] * ":", organism)
  end
end

Workflow.require_workflow "Appris"

require 'rbbt/entity/appris'

post '/wizard' do
  template_render('wizard', params, "Wizard", :cache_type => :sync)
end

get '/get_pdb' do
  query = consume_parameter :query
  content_type 'text/plain'
  send_file Open.open(query)
end

get '/appris_features' do
  isoform = consume_parameter :isoform
  raise "No isoform provided" if isoform.nil?
  features = Proteomics.appris_features(isoform)
  content_type :json
  halt 200, features.to_json
end

get '/interactome3d' do
  puts '/interactome3d'
  pdb = consume_parameter :pdb
  type1 = consume_parameter :type1
  type2 = consume_parameter :type2
  url ="http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=#{type1}&type2=#{type2}&pdb=#{pdb}"
  raise "No interactome3d provided" unless [pdb, type1, type2].all?
   content_type 'text/plain'
  send_file Open.open(url)
end
