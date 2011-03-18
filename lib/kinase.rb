require 'rbbt-util'
require 'rbbt/util/workflow'
require 'rbbt/util/cmd'
module Kinase
  extend WorkFlow

  ROOT = File.dirname(File.dirname(File.expand_path(__FILE__)))
  BIN_DIR = File.join(ROOT, 'bin')

  @basedir = File.join(ROOT, 'var')

  task_option :list, "Lista de mutaciones", :string
  task :input => :string do |list|
    list
  end

  task :patterns => :string do
    result = CMD.cmd("perl -I #{BIN_DIR} #{File.join(BIN_DIR, 'PatternGenerator.pl')} #{ previous_jobs["input"].path } #{File.join(ROOT, "etc/feature.number.list")}").read
    result
  end

  task :predict => :string do 
    CMD.cmd(File.join(BIN_DIR, "run_svm.py --m=e --o=#{File.join(task.workflow.basedir, task.name)} \
    --svm=#{File.join(ROOT, 'share/model/final.svm')} --cfg=#{File.join(ROOT, 'etc/svm.config')}"), 
    "--ts=" => previous_jobs["patterns"].path)

    nil
  end
end

#puts Kinase.job(:predict, "test", Open.read(File.join(Kinase::ROOT, 'data/EXAMPLES/test.input'))).run.load
