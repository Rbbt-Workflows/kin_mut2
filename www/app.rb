require File.join(Sinatra::Application.root,'../lib/kinase')
require 'rbbt/util/tsv'
require 'pp'


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
    @filtered = job.input("patterns").info[:filtered_out]

    haml :result
  end
end

get '/filtered/:name' do
  job = Kinase.load_job(:predict, params[:name])

  content_type "text/plain"
  job.input("patterns").info[:filtered_out]
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
