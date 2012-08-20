require File.join(Sinatra::Application.root,'../lib/kinase')
require 'rbbt/tsv'
require 'rbbt/sources/entrez'
require 'rbbt/sources/go'
require 'rbbt/sources/pfam'
require 'rbbt/sources/uniprot'
require 'rbbt/entity/gene'
require 'pp'

$kinase_groups = Kinase.data["kinase_group_description.tsv"].find(:lib).tsv :single
$prot_goterms = Kinase.data["uniprot2go.txt"].find(:lib).tsv :single, :fields => [2]

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
  @title = "wKinMut: Kinase Mutation Damange Prediction"
  haml :help
end

get '/job/:name' do
  @title = "wKinMut: #{params[:name].match(/(.*)_(.*)/)[1] }"
  job = Kinase.load_id(File.join("default", params["name"]))

  while not job.done?
    job.join
    job = Kinase.load_id(File.join("default", params["name"]))
  end

  if job.error?
    @message = "Error in job #{ job.name }: #{job.messages.last}"
    haml :error
  else
    begin
      @res = TSV.open(job.step(:predict).path, :key_field => 2, :sep => /\s+/)
      @res.unnamed = true
      @job = params[:name]
      @jobname, @hash = params[:name].match(/(.*)_(.*)/).values_at 1, 2
      @translations = job.step("patterns").step("input").info[:translations]
      @translations_id = job.step("patterns").step("input").info[:translations_id]
      @list = job.step("input").info[:inputs][:list].split("\n")
      @filtered = (job.step("patterns").info[:filtered_out] || []) + job.step("patterns").step("input").info[:synonymous]
      @uniprot_groups = {}

      @res.each{|mutation,values|
        mutation = mutation[1..-1]
        prot, m = mutation.split(/_/)
        @uniprot_groups[mutation] = Kinase.get_features(job, prot, m)["uniprot_group"]
      }

      haml :result
    rescue
      if $!.message == "Empty content"
        @message = "Error in job #{ job.name }: No results produced, maybe no kinases where identified in the input."
        haml :error
      end
      raise $!
    end
  end
end

get '/filtered/:name' do
  job = Kinase.load_id(File.join("default", params[:name]))

  @translations = job.step("input").info[:translations]

  @filtered = job.step("patterns").info[:filtered_out]
  @synonymous = job.step("input").info[:synonymous] 

  haml :missing
end


get '/original/:name' do
  job = Kinase.load_id(File.join("default", params[:name]))

  content_type "text/plain"
  job.step("input").info[:inputs][:list]
end

get '/download/:name' do
  job = Kinase.load_id(File.join("default", params[:name]))

  content_type "text/plain"
  res = TSV.open job.step(:predict).path, :key_field => 2, :sep => /\s+/

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
  job = Kinase.load_id(File.join("default", params[:name]))
  @protein, @mutation = params.values_at :protein, :mutation

  @title = "wKinMut: #{params[:name].match(/(.*)_(.*)/)[1] } > #{@protein} > #{@mutation}"

  entrez_index = Organism::Hsa.identifiers.index(:target => "Entrez Gene ID", :persist =>  true)
  name_index   = Organism::Hsa.identifiers.index(:target => "Associated Gene Name", :persist =>  true)
  ensp_index   = Organism::Hsa.protein_identifiers.index(:target => "Ensembl Protein ID", :persist =>  true)
  ensg_index   = Organism::Hsa.identifiers.index(:target => "Ensembl Gene ID", :persist =>  true)

  @translations = job.step("patterns").step("input").info[:translations]
  @translations_id = job.step("patterns").step("input").info[:translations_id]

  @res = TSV.open(job.step(:predict).path, :key_field => 2, :sep => /\s+/)

  @ensp = (ensp_index[@protein] || []).first
  @ensg = (ensg_index[@protein] || []).first
  @name = (name_index[@ensp] || []).first
  @entrez = (entrez_index[@name] || []).first
  @jobname = params[:name]
  @other = job.step(:other_predictors).load

  unless @entrez.nil?
    gene = Entrez.get_gene(@entrez)
    @description = gene.description
    @summary     = gene.summary
  end

  @goterms = Misc.zip_fields(Organism.gene_go("Hsa").tsv(:persist => true, :unnamed => true)[@ensg])
  @features = Kinase.get_features(job, @protein, @mutation)

  haml :details
end


post '/' do

  if params[:file] && params[:file][:tempfile]
    mutations = params[:file][:tempfile].read
  else
    mutations = params[:mutations].upcase
  end

  jobname = params[:jobname].gsub(/\s+/,'_')
  jobname = "JOB" if jobname.nil? or jobname.empty?
  job = Kinase.job(:default, jobname , :list =>  mutations)
  job.fork
  if params[:jobid]
    job.name
  else
    redirect "/job/#{job.name}"
  end
end

get '/' do
  haml :index
end

get '/jmol/:uniprot/:position' do
  uniprot = params[:uniprot]
  position = params[:position]

  pdbs = Kinase.pdb_position(uniprot, position)
  haml :jmol, :layout => false, :locals => {:pdbs => pdbs}
end

get '/sentences/:uniprot' do
  uniprot = params[:uniprot]

  haml :sentences, :layout => false, :locals => {:uniprot => uniprot}
end
