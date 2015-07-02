require 'kinase'
require 'rbbt/tsv'
require 'rbbt/sources/entrez'
require 'rbbt/sources/go'
require 'rbbt/sources/pfam'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/string'
require 'rbbt/entity/gene'
require 'rbbt/entity/gene/basic'
require 'rbbt/entity/protein'
require 'rbbt/entity/protein/basic'
require 'pp'

Workflow.require_workflow "Genomics"
Workflow.require_workflow "Structure"
Workflow.require_workflow "KinaseSARfari"

require 'rbbt/entity/mutated_isoform'


add_workflow Structure, true

require 'structure/uniprot'

organism = KinMut2.organism

$protein_kinase_groups = Rbbt.data["Uniprot.kinase.basicinfo.txt"].tsv(:key_field => "UNI_ACC", :fields => ["GROUP"], :persist => true, :type => :single)
$kinase_groups = Kinase.data["kinase_group_description.tsv"].find(:lib).tsv :single
$kinase_groups_lor = Rbbt.data["Uniprot.kinase.groups_lor.txt"].tsv(:key_field => "group", :fields => ["logodds"], :persist => true, :type => :single, :header_hash => "")
$GO_lor = Rbbt.data["Gene_Ontology.GOLOR.txt"].tsv(:key_field => "GOterm", :fields => ["golor"], :persist => true, :type => :single, :header_hash => "")

$prot_goterms = Kinase.data["uniprot2go.txt"].find(:lib).tsv :single, :fields => [2]
$domains_lor = Rbbt.data["Uniprot.kinase.domains_lor.txt"].tsv(:key_field => "domain", :fields => ["logodds"], :persist => true, :type => :single, :header_hash => "")
$gene_essentiality = Rbbt.data["gene_essentiality"].tsv :type => :single

$goterm_score = Kinase.data["GOlogoddsratio.perterm.txt"].find(:lib).tsv :single, :fields => ["DESCRIPTION" ]
$go_terms = Organism.gene_go(organism).tsv(:persist => true, :unnamed => true)

$organism = Organism.default_code(organism)

$pfam_names = {}

$pfam_info = Pfam.domains.tsv 
$pfam_nam2ac = {}

ac = nil
description = ""
content = ""
Open.read(Kinase.data["Pfam.desc"].find) do |line|
  code, value = line.chomp.match(/^(\w\w)\s+(.*)/).values_at 1,2
  case code
  when "AC"
    $pfam_names[ac] = {:description => description, :content => content} unless ac.nil?
    $pfam_nam2ac[description] = ac
    ac = value.sub(/\..*/,'')
    description = ""
    content = ""
  when "DE"
    description += value
  when "CC"
    content += value
  end
end


get '/job/:name' do
  @title = "wKinMut: #{params[:name].match(/(.*)_(.*)/)[1] }"
  job = Kinase.load_id(File.join("default", params["name"]))

  sleep 3 unless job.done?

  while not job.done? and not job.error?
    return template_render("wait", :job => job.name)
  end

  if job.error?
    if job.messages[-2] =~ /No predictions/
      @job = job.step(:patterns)
      template_render("no_predicitons")
    else
      @message = "Error in job #{ job.name }: #{job.messages[-2]}"
      template_render("error")
    end
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

      template_render("result")
    rescue
      if $!.message == "Empty content"
        @message = "Error in job #{ job.name }: No results produced, maybe no kinases where identified in the input."
        template_render("error")
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
  @job = KinMut2.load_id(File.join("predict_all", params[:name]))
  @protein, @mutation = params.values_at :protein, :mutation
  @position = @mutation.scan(/\d+/).first.to_i

  @title = "wKinMut: #{params[:name].match(/(.*)_(.*)/)[1] } > #{@protein} > #{@mutation}"

  @ensp = $ensp_index[@protein] || $ensp_index_all[@protein] || [].first
  @ensg = $ensg_index[@protein] || [].first
  @name = $name_index[@protein] || [].first
  @entrez = $entrez_index[@protein] || [].first

  if @ensp and @ensp.sequence
    @map = Structure.uniprot_sequence_map(@protein, @ensp.sequence) 
    @fixed_position  = @map[@position]
    @fixed_mutation  = @mutation.sub(@position.to_s, @fixed_position.to_s)
  else
    @map = nil
    @fixed_position  = @position
  end

  @jobname = params[:name]

  unless @entrez.nil?
    gene = Entrez.get_gene(@entrez)
    @description = gene.description
    @summary     = gene.summary
  end

  @goterms = Misc.zip_fields($go_terms[@ensg])

  template_render("details", {}, "Details #{[@job.name, @protein, @mutation] * "-"}", :cache_type => :async)
end

post '/job' do
  if params[:file] && params[:file][:tempfile]
    mutations = params[:file][:tempfile].read
  else
    mutations = params[:mutations].upcase
  end

  mutations = mutations.gsub(/\b +\b/,'_')
  
  puts mutations

  jobname = params[:jobname].gsub(/\s+/,'_')
  jobname = "JOB" if jobname.nil? or jobname.empty?
  job = Kinase.job(:default, jobname , :list =>  mutations)
  job.fork
  sleep 1
 
  halt 200, job.name
end

get '/status/:name' do
  job = Kinase.load_id(File.join("default", params[:name]))
  halt 200, job.status.to_s
end

post '/' do

  if params[:file] && params[:file][:tempfile]
    mutations = params[:file][:tempfile].read
  else
    mutations = params[:mutations].upcase
  end

  mutations = mutations.gsub(/\b +\b/,'_')
  
  puts mutations

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

get '/input' do
  template_render("input")
end

get '/resources' do
  template_render("resources")
end

get '/references' do
  template_render("references")
end

get '/jmol/:uniprot/:position' do
  uniprot = params[:uniprot]
  position = params[:position]

  pdbs = Kinase.pdb_position(uniprot, position)
  template_render("jmol", :_layout => false, :pdbs => pdbs)
end

get '/sentences/:uniprot' do
  uniprot = params[:uniprot]

  template_render("sentences", :uniprot => uniprot)
end
