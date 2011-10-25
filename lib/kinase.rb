require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/resource'
require 'rbbt/util/cmd'
require 'rbbt/sources/organism'
require 'nokogiri'
require 'pg'


Workflow.require_workflow 'translation'
module Kinase
  extend Workflow
  extend Resource

  class << self
    include LocalPersist
    local_persist_dir = File.join(File.dirname(__FILE__), '../cache')
  end


  Kinase.software.opt.svm_light.claim :install, Rbbt.share.install.software.svm_light.find

  module Postgres
    HOST = 'padme'
    PORT = nil
    OPTIONS = nil
    TTY = nil
    DBNAME = 'tm_kinase_muts'

    def self.driver
      PGconn.connect(:host => HOST, :port => PORT, :dbname => DBNAME, :user => 'jmgonzalez')
    end

    def self.snp2l(uniprot, position)
      query =<<-EOT
SELECT
tum.wt_aa,
tum.mutant_aa,
tum.pubmed_id,
tum.new_position,
tum.old_position,
td.type,
mm.type,
ts.line_content

FROM 
tm_updated_mentions tum,
tm_datasets td,
mapping_methods mm,
tm_sentences ts

WHERE 
tum.tm_dataset_id=td.id AND
tum.mapping_method_id=mm.id AND
tum.acc='#{uniprot}' AND
tum.new_position=#{position} AND
tum.oldtriplet=ts.oldtriplets

ORDER BY 
tum.new_position,
tum.wt_aa, 
tum.mutant_aa
;
      EOT
      res = driver.exec(query)
      res
    end

    def self.snp2db(uniprot, position)
      query =<<-EOT
SELECT 
cm.acc,
cm.seq_pos,
cm.wt,
cm.mutant,
cm.description,
ed.type 
FROM 
complementary_muts cm,
external_dbs ed
WHERE 
cm.external_db_id = ed.id AND
acc='#{uniprot}' AND
seq_pos=#{position}
;
      EOT
      res = driver.exec(query)
      res
    end
  end

  def self.ihop_interactions(uniprot)
    url = "http://ws.bioinfo.cnio.es/iHOP/cgi-bin/getSymbolInteractions?ncbiTaxId=9606&reference=#{uniprot}&namespace=UNIPROT__AC" 
    doc = Nokogiri::XML(Open.read(url))
    sentences = doc.css("iHOPsentence")
    sentences
  end

  def self.error_in_wt_aa?(protein, mutation)
    wt, pos, m = mutation.match(/([A-Z])(\d+)([A-Z])/i).values_at 1,2,3

    @@sequences ||= self.local_persist("sequence", :tsv, :source => data["KinaseAccessions_Group_Seqs.txt"].find) do 
      TSV.open data["KinaseAccessions_Group_Seqs.txt"].find, :type => :single, :fields => [2]
    end

    real_wt = @@sequences[protein][pos.to_i - 1].chr

    if wt == real_wt
      false
    else
      real_wt
    end
  end

  def self.get_features(job, protein, mutation)
    @@feature_names ||= self.etc["feature.number.list"].find(:lib).tsv(:type => :single).sort_by{|key,value| key.to_i}.collect{|key, value| value}

    @@patterns ||= {}

    #patterns = @@patterns[job] ||= TSV.open(job.step("patterns").path, :list, :fields => @@feature_names, :key_field => @@feature_names.length, :fix => Proc.new{|l| l.sub('#','').sub(/^\d+\t/,'').gsub(/\d+:/,'')})
    patterns = @@patterns[job] ||= TSV.open(job.step("patterns").path, :list, :key_field => @@feature_names.length, :fix => Proc.new{|l| l.sub('#','').sub(/^\d+\t/,'').gsub(/\d+:/,'')})
    patterns.key_field = "Protein Mutation"
    patterns.fields = @@feature_names

    pattern = patterns[[protein, mutation] * "_"]
    info = {}

    %w(SIFTscore SIFTscore_binned TDs_fscore_diff TDs_fscore_mt TDs_fscore_wt
    biochem_diffkdhydrophobicity firedb pfam_any phosphoelm sumGOLOR
    swannot_act_site swannot_act_site swannot_any swannot_binding
    swannot_carbohyd swannot_catalytic swannot_disulfid swannot_metal
    swannot_mod_res swannot_mutagen swannot_np_bind swannot_ptm swannot_signal
    swannot_site swannot_transmem).each do |key|
      info[key] = pattern[key]
    end

    info["uniprot_group"] = %w( class_uniprotgroup_AGC class_uniprotgroup_Atypical_ADCK
    class_uniprotgroup_Atypical_Alpha-type class_uniprotgroup_Atypical_FAST
    class_uniprotgroup_Atypical_PDK-BCKDK class_uniprotgroup_Atypical_PI3-PI4
    class_uniprotgroup_Atypical_RIO class_uniprotgroup_CAMK
    class_uniprotgroup_CK1 class_uniprotgroup_CMGC class_uniprotgroup_NEK
    class_uniprotgroup_Other class_uniprotgroup_RGC class_uniprotgroup_STE
    class_uniprotgroup_TK class_uniprotgroup_TKL).select{|key|
      pattern[key] == "1"
    }.first
    info["uniprot_group"].sub!(/class_uniprotgroup_/,'') unless info["uniprot_group"].nil?

    info["pfam"] = %w( pfam_PF00017 pfam_PF00018 pfam_PF00023 pfam_PF00027 pfam_PF00028
    pfam_PF00041 pfam_PF00047 pfam_PF00051 pfam_PF00063 pfam_PF00130
    pfam_PF00168 pfam_PF00169 pfam_PF00211 pfam_PF00226 pfam_PF00412
    pfam_PF00415 pfam_PF00433 pfam_PF00435 pfam_PF00454 pfam_PF00498
    pfam_PF00520 pfam_PF00531 pfam_PF00536 pfam_PF00560 pfam_PF00564
    pfam_PF00566 pfam_PF00567 pfam_PF00595 pfam_PF00611 pfam_PF00612
    pfam_PF00615 pfam_PF00621 pfam_PF00629 pfam_PF00659 pfam_PF00754
    pfam_PF00757 pfam_PF00779 pfam_PF00780 pfam_PF00787 pfam_PF01030
    pfam_PF01064 pfam_PF01094 pfam_PF01163 pfam_PF01392 pfam_PF01403
    pfam_PF01404 pfam_PF01833 pfam_PF02019 pfam_PF02185 pfam_PF02259
    pfam_PF02816 pfam_PF02828 pfam_PF03109 pfam_PF03607 pfam_PF03623
    pfam_PF06293 pfam_PF06479 pfam_PF07647 pfam_PF07686 pfam_PF07699
    pfam_PF07701 pfam_PF07714 pfam_PF08064 pfam_PF08163 pfam_PF08238
    pfam_PF08311 pfam_PF08332 pfam_PF08368 pfam_PF08477 pfam_PF08515
    pfam_PF08919 pfam_PF08926 pfam_PF09027 pfam_PF09042 pfam_PF10409
    pfam_PF10436 pfam_PF11555 pfam_PF11640 pfam_PF12063 pfam_PF12179
    pfam_PF12202 pfam_PF12474).select{|key|
      pattern[key] == "1"
    }.collect{|v| v.sub(/pfam_/,'')} * "|"

    info
  end

  input :list, :string, "Lista de mutations"
  task :input => :string do |list|
    proteins = []
    mutations = []

    list.split(/\n/).each{|l| 
      l.strip!
      if l.match(/(.*)[_ \t,]+(.*)/)
        prot, mut = $1, $2
        proteins << prot
        mutations << mut
      end
    }

    set_info :originals, proteins

    translated = Translation.job(:translate_protein, "", :proteins => proteins, :format => "UniProt/SwissProt Accession").exec

    set_info :translated, translated

    translations = Hash[*(translated.zip(proteins)).flatten]
    
    set_info :translations, translations

    list = translated.zip(mutations)

    same_aa = list.select{|p,m| m[0] == m[-1]}

    set_info :synonymous, same_aa

    list.reject!{|p,m| m[0] == m[-1]}

    list.reject{|p,m| p.nil?}.collect{|p,m| [p,m] * "_"} * "\n"
  end

  dep :input
  task :patterns => :string do
    error_file = TmpFile.tmp_file
    patterns = CMD.cmd("perl -I #{Kinase.bin.find} #{Kinase['bin/PatternGenerator.pl'].find} #{ step("input").path } #{Kinase["etc/feature.number.list"].find} 2> #{error_file}").read
    if Open.read(error_file).any? and Open.read(error_file) =~ /is not a valid/
        set_info :filtered_out, Open.read(error_file).split(/\n/).collect{|l| l.match(/(\w*) is not a valid/)[1]}
    else
      set_info :filtered_out, []
    end

    patterns
  end

  dep :patterns
  task :predict => :string do 
    CMD.cmd("#{Kinase["bin/run_svm.py"].find(:lib)} --m=e --o=#{path}.files \
    --svm=#{Kinase['share/model/final.svm'].find} --cfg=#{Kinase['etc/svm.config'].find}", 
    "--ts=" => step("patterns").path)
    FileUtils.mv File.join(path + '.files', File.basename(step("patterns").path)), path

    nil
  end
end

#puts Kinase.job(:predict, "test", Open.read(File.join(Kinase::ROOT, 'data/EXAMPLES/test.input'))).run.load

if __FILE__ == $0

  puts Kinase.ihop_interactions("P07949").first

  exit
  ddd Kinase::Postgres.snp2l('P07949', 806).to_a

  exit
  job = Kinase.job(:predict, "test", :list => Open.read(File.join(['data/EXAMPLES/test.input'])))
  job.clean.run
  ddd job

  ddd Kinase.get_features(job, 'O14936', 'P396S')

  
end
