require 'rbbt-util'
require 'rbbt/util/workflow'
require 'rbbt/util/cmd'
require 'rbbt/sources/organism'

module Kinase
  extend WorkFlow

  Rbbt.software.opt.svm_light.define_as_install Rbbt.share.install.software.svm_light.find

  def self.error_in_wt_aa?(protein, mutation)
    wt, pos, m = mutation.match(/([A-Z])(\d+)([A-Z])/).values_at 1,2,3

    sequences = local_persist(data["KinaseAccessions_Group_Seqs.txt"], :TSV, :tsv) do |file, *other|
      TSV.new Open.open(file), :single, :fields => 2
    end

    real_wt = sequences[protein][pos.to_i - 1].chr

    if wt == real_wt
      false
    else
      real_wt
    end
  end

  task_option :list, "Lista de mutations", :string
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

    translated = Organism::Hsa.normalize(proteins, "UniProt/SwissProt Accession")

    set_info :translated, translated

    translations = Hash[*(translated.zip(proteins)).flatten]
    
    set_info :translations, translations

    list = translated.zip(mutations)

    ddd list
    same_aa = list.select{|p,m| m[0] == m[-1]}

    set_info :synonymous, same_aa

    list.reject!{|p,m| m[0] == m[-1]}

    list.reject{|p,m| p.nil?}.collect{|p,m| [p,m] * "_"} * "\n"
  end

  task :patterns => :string do
    error_file = TmpFile.tmp_file
    patterns = CMD.cmd("perl -I #{Kinase.bin.find} #{Kinase['bin/PatternGenerator.pl'].find} #{ previous_jobs["input"].path } #{Kinase["etc/feature.number.list"].find} 2> #{error_file}").read
    if Open.read(error_file).any?
      set_info :filtered_out, Open.read(error_file).split(/\n/).collect{|l| l.match(/(\w*) is not a valid/)[1]}
    else
      set_info :filtered_out, []
    end

    patterns
  end

  task :predict => :string do 
    CMD.cmd("#{Kinase["bin/run_svm.py"].find} --m=e --o=#{File.join(Kinase.jobdir, task.name)} \
    --svm=#{Kinase['share/model/final.svm'].find} --cfg=#{Kinase['etc/svm.config'].find}", 
    "--ts=" => previous_jobs["patterns"].path)
    FileUtils.mv File.join(Kinase.jobdir, task.name, File.basename(previous_jobs["patterns"].path)), path

    nil
  end
end

#puts Kinase.job(:predict, "test", Open.read(File.join(Kinase::ROOT, 'data/EXAMPLES/test.input'))).run.load

if __FILE__ == $0

  puts Kinase.check_wt_aa('P05129', 'A523D')

end
