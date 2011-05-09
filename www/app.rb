require File.join(Sinatra::Application.root,'../lib/kinase')
require 'rbbt/util/tsv'
require 'rbbt/sources/entrez'
require 'rbbt/sources/go'
require 'pp'

$kinase_groups = Kinase.data["kinase_group_description.tsv"].tsv :single
$prot_goterms = Kinase.data["uniprot2go.txt"].tsv :single, :fields => 2
$goterm_score = Kinase.data["GOlogoddsratio.perterm.txt"].tsv :single, :fields => ["DESCRIPTION" ]

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
  haml :help
end

get '/job/:name' do
  job = Kinase.load_job(:predict, params[:name])
  if not job.done?
    job.join
    job.load_dependencies
  end

  if job.error?
    "Error in job #{ job.id }: #{job.messages.last}"
  else
    @res = TSV.new job.open, :key => 2, :sep => /\s+/

    @job = params[:name]
    @jobname, @hash = params[:name].split(/_/)
    @translations = job.input("input").info[:translations]
    @filtered = (job.input("patterns").info[:filtered_out] || []) + job.input("input").info[:synonymous]
    @uniprot_groups = {}
    @res.each{|mutation,values|
      mutation = mutation.sub('#', '')
      prot, m = mutation.split(/_/)
      @uniprot_groups[mutation] = Kinase.get_features(job, prot, m)["uniprot_group"]
    }

    haml :result
  end
end

get '/filtered/:name' do
  job = Kinase.load_job(:predict, params[:name])

  @translations = job.input("input").info[:translations]
  
  @filtered = job.input("patterns").info[:filtered_out]
  @synonymous = job.input("input").info[:synonymous] 
  
  haml :missing
end


get '/original/:name' do
  job = Kinase.load_job(:predict, params[:name])

  content_type "text/plain"
  job.input("input").info[:options][:list]
end

get '/download/:name' do
  job = Kinase.load_job(:predict, params[:name])

  content_type "text/plain"
  res = TSV.new job.open, :key => 2, :sep => /\s+/
  translations = job.input("input").info[:translations]

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
  job = Kinase.load_job(:predict, params[:name])
  @protein, @mutation = params.values_at :protein, :mutation

  index = Organism::Hsa.identifiers.index(:target => "Entrez Gene ID", :persistence =>  true)
  name = Organism::Hsa.identifiers.index(:target => "Associated Gene Name", :persistence =>  true)

  @entrez = (index[@protein] || []).first
  @name = (name[@protein] || []).first

  unless @entrez.nil?
    @description = Entrez.get_gene(@entrez).description
    @summary     = Entrez.get_gene(@entrez).summary
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
  job = Kinase.job(:predict, jobname , mutations)
  job.fork
  redirect "/job/#{job.id}"
end

get '/' do
  haml :index
end
