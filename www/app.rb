require '../lib/kinase'
require 'rbbt/util/tsv'

get '/job/:name' do
  job = Kinase.job(:predict, params[:name], Open.read(File.join(Kinase::ROOT, 'data/EXAMPLES/test.input').sub(/\n.*/ms,'\n'))).fork
  job.join
  puts job.read
  @res = TSV.new job.open, :key => 2, :sep => /\s+/
    puts @res.to_s
  
  haml :result
end

post '/' do
  session[:list] = 
  session[:job] = job = Kinase.job(:predict, "test", Open.read(File.join(Kinase::ROOT, 'data/EXAMPLES/test.input')))
  job.fork
  redirect "/job/#{job.name}"
end

get '/' do

end
