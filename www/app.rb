require File.join(Sinatra::Application.root,'../lib/kinase')
require 'rbbt/tsv'
require 'rbbt/sources/entrez'
require 'rbbt/sources/go'
require 'pp'

$kinase_groups = Kinase.data["kinase_group_description.tsv"].find(:lib).tsv :single
$prot_goterms = Kinase.data["uniprot2go.txt"].find(:lib).tsv :single, :fields => 2
$goterm_score = Kinase.data["GOlogoddsratio.perterm.txt"].find(:lib).tsv :single, :fields => ["DESCRIPTION" ]

$pfam_names = {}

ac = nil
description = ""
content = ""
Open.read(Kinase.data["Pfam.desc"].find) do |line|
  code, value = line.chomp.match(/^(\w\w)\s+(.*)/).values_at 1,2
  case code
  when "AC"
    $pfam_names[ac] = {:description => description, :content => content} unless ac.nil?
    ac = value.sub(/\..*/,'')
    description = ""
    content = ""
  when "DE"
    description += value
  when "CC"
    content += value
  end
end

get '/help' do
  @title = "KinMut: Kinase Mutation Damange Prediction"
  haml :help
end

get '/job/:name' do
  @title = "KinMut: #{params[:name].sub(/_.*/,'') }"
  job = Kinase.load_id(File.join("predict", params["name"]))

  while not job.done?
    job.join
    job = Kinase.load_id(File.join("predict", params["name"]))
  end

  if job.error?
    "Error in job #{ job.name }: #{job.messages.last}"
  else
    @res = TSV.open job.path, :key_field => 2, :sep => /\s+/

    @job = params[:name]
    @jobname, @hash = params[:name].split(/_/)
    @translations = job.step("patterns").step("input").info[:translations]
    @list = job.step("input").info[:inputs][:list].split("\n")
    @filtered = (job.step("patterns").info[:filtered_out] || []) + job.step("patterns").step("input").info[:synonymous]
    @uniprot_groups = {}
    @res.each{|mutation,values|
      mutation = mutation.sub('#', '')
      prot, m = mutation.split(/_/)
      @uniprot_groups[mutation] = Kinase.get_features(job, prot, m)["uniprot_group"]
    }

    index = Organism::Hsa.identifiers.index(:target => "Entrez Gene ID", :persist =>  true)

    Entrez.get_gene(@res.keys.collect{|mutation|  
      prot, m = mutation.sub('#','').split(/_/)
      (index[prot] || []).first
    }.compact)

    haml :result
  end
end

get '/filtered/:name' do
  job = Kinase.load_id(File.join("predict", params[:name]))

  @translations = job.step("input").info[:translations]
  
  @filtered = job.step("patterns").info[:filtered_out]
  @synonymous = job.step("input").info[:synonymous] 
  
  haml :missing
end


get '/original/:name' do
  job = Kinase.load_id(File.join("predict", params[:name]))

  content_type "text/plain"
  job.step("input").info[:options][:list]
end

get '/download/:name' do
  job = Kinase.load_id(File.join("predict", params[:name]))

  content_type "text/plain"
  #res = TSV.new job.open, :key => 2, :sep => /\s+/
  res = TSV.open job.path, :key_field => 2, :sep => /\s+/

  translations = job.step("input").info[:translations]

  line = []
  res.collect do |mutation,values|
    line = []
    mutation = mutation.sub('#', '')
    prot, m = mutation.split(/_/)
    if translations[prot] != prot
      line << translations[prot] << " " << prot
    else
      line << prot
    end
    
    line << m

    line << values[1].first.to_f
    if values[1].first.to_f < -0.5  
      line << "Neutral"
    else
      line << "Damaging"
    end
    
    line * "\t"
  end * "\n"
end

get '/details/:name/:protein/:mutation' do
  job = Kinase.load_id(File.join("predict", params[:name]))
  @protein, @mutation = params.values_at :protein, :mutation

  @title = "KinMut: #{params[:name].sub(/_.*/,'') } > #{@protein} > #{@mutation}"

  index = Organism::Hsa.identifiers.index(:target => "Entrez Gene ID", :persist =>  true)
  name = Organism::Hsa.identifiers.index(:target => "Associated Gene Name", :persist =>  true)

  @entrez = (index[@protein] || []).first
  @name = (name[@protein] || []).first

  unless @entrez.nil?
    gene = Entrez.get_gene(@entrez)
    @description = gene.description
    @summary     = gene.summary
  end

  @goterms = Misc.process_to_hash $prot_goterms[@protein].split(/;/) do |list| list.collect{|id| $goterm_score[id]} end

  @features = Kinase.get_features(job, @protein, @mutation)

  haml :details
end


post '/' do

  if params[:file] && params[:file][:tempfile]
    mutations = params[:file][:tempfile].read
  else
    mutations = params[:mutations]
  end

  jobname = params[:jobname]
  jobname = "JOB" if jobname.nil? or jobname.empty?
  job = Kinase.job(:predict, jobname , :list =>  mutations)
  job.fork
  redirect "/job/#{job.name}"
end

get '/' do
  haml :index
end

get '/sentences/:uniprot' do
  uniprot = params[:uniprot]

  haml :sentences, :layout => false, :locals => {:uniprot => uniprot}
end
