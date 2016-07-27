require 'rbbt/workflow'
require 'rbbt/util/misc/exceptions'
require 'rbbt/sources/clinvar'

Workflow.require_workflow "Genomics"
Workflow.require_workflow "Structure"
Workflow.require_workflow "DbNSFP"


module KinMut2
  extend Workflow

  helper :predict_script do
    Rbbt.bin["predict_mutations.pl"].find(:lib)
  end

  def self.organism
    Organism.default_code("Hsa")
  end

  def self.all_kinases
    Kinase.data["Uniprot.kinase.accessions.txt"].find(:lib).read.split("\n")
  end


  $entrez_index = Organism.identifiers(organism).index(:target => "Entrez Gene ID", :fields => ["UniProt/SwissProt Accession"], :persist =>  true, :order => true)
  $name_index   = Organism.identifiers(organism).index(:target => "Associated Gene Name", :fields => ["UniProt/SwissProt Accession"], :persist =>  true, :order => true)
  $ensg_index   = Organism.identifiers(organism).index(:target => "Ensembl Gene ID", :fields => ["UniProt/SwissProt Accession"], :persist =>  true, :order => true)
  $ensp_index_all = Organism.protein_identifiers(organism).index(:target => "Ensembl Protein ID", :fields => ["UniProt/SwissProt Accession"], :persist =>  true,  :order => true).tap{|i| i.namespace = organism}
  $ensp_index   = Organism.protein_identifiers(organism).index(:target => "Ensembl Protein ID", :fields => ["UniProt/SwissProt Accession"], :persist =>  true,  :order => true, :data_tsv_grep => Appris::PRINCIPAL_ISOFORMS.to_a)
  $kinase_FDA_drugs = Rbbt.data["FDA_drugs_kinases_all_available.txt"].tsv :persist => false, :header_hash => "", :fields => %w(Sponsor Indications Target), :sep2 => /\s*,\s+/
  $clinvar = ClinVar.mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true

  helper :organism do
    KinMut2.organism
  end


  input :mutations, :array, "Mutations"
  task :translate => :array do |mutations|
    raise ParameterException, "No mutations provided" if mutations.nil? 

    translated = []
    organism = Organism.default_code("Hsa")
    translations = TSV.setup({}, :key_field => "Protein", :fields => ["UniProt/SwissProt Accession"], :type => :single, :namespace => organism)

    uni = Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :fields => ["Ensembl Protein ID"], :persist => true)
    name = Organism.identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true)
    missing = []
    skipped_mutations = []
    mutation_translations = TSV.setup({}, :key_field => "Final", :fields => ["Orginal"], :type => :single)
    TSV.traverse mutations, :into => translated do |mutation|
      mutation = mutation.first if Array === mutation
      protein, change = mutation.split(/[,\s:]+/)
      if change =~ /^([A-Z])\d+([A-Z])$/ and not $1 == $2
        translation = case protein
                      when /ENSP/
                        translations[protein] = uni[protein]
                      else
                        translations[protein] = name[protein] || protein
                      end
        if translation.nil?
          missing << protein
          skipped_mutations << mutation
          next
        end
        mutation_translations[[translation, change] * " "] = mutation
        [translation, change] * " "
      else
        skipped_mutations << mutation
        next
      end
    end

    Open.write(file(:translations), translations.to_s)
    Open.write(file(:mutation_translations), mutation_translations.to_s)
    Open.write(file(:skipped_mutations), skipped_mutations*"\n")
    set_info :missing_translations, missing

    translated
  end

  dep :translate
  task :predict => :tsv do 
    mutations = step(:translate).load
    script = self.predict_script

    output = files_dir
    TmpFile.with_file(mutations*"\n") do |input|
      `perl '#{script}' -input '#{input}' -output '#{output}'`
    end
  
    features_tsv = {}
    arff = file('vectors.weka.arff').read
    attributes = arff.scan(/@attribute (.*) .*/).flatten
    vectors = begin
                arff.split('@data').last.split("\n")[1..-1].collect{|l| l.split ","}
              rescue
                raise ParameterException, "No valid kinase mutations" 
              end
    mutations.zip(vectors).each do |mutation,vector|
      features_tsv[mutation.sub(/\s/,' ')] = vector
    end
    TSV.setup(features_tsv, :key_field => "Mutation", :fields => attributes, :type => :list)

    Open.write(file("features"), features_tsv.to_s)

    FileUtils.cp(step(:translate).file(:translations), file(:translations))

    fields = %w(Mutation Prediction Score)
    header = "#" << fields * "\t"
    file_txt = header << "\n" << Open.read(File.join(output, "vectors.weka.predictions"))
    tsv = TSV.open(file_txt)
    #index = TSV.index(file(:translations), :target => "Protein", :merge => true)
    index = file(:translations).tsv :fields => ["Protein"],:type => :flat, :merge => true, :key_field => "UniProt/SwissProt Accession"
    tsv = tsv.to_double.add_field "Fixed mutation" do |mutation,values|
      uni, change = mutation.split(" ")
      origs = index[uni]
      origs.collect do |orig|
        [orig, change]  * ":"
      end
    end

    Open.write(file(:pre_fixed), tsv.to_s)
    
    tsv = tsv.reorder("Fixed mutation", ["Prediction", "Score"])
    tsv.key_field = "Mutated Isoform"
    Open.write(file(:fixed), tsv.to_s)

    file_txt
  end

  dep :predict
  task :predict_fix => :tsv do
    TSV.get_stream step(:predict).file(:fixed)
  end

  dep :predict
  task :predict_ensp => :tsv do
    tsv = step(:predict).load
    fields = tsv.fields
    tsv.key_field = "Original Mutated Isoform"
    tsv.add_field "Mutated Isoform" do |mi,values|
      prot, change = mi.split(" ")
      ensp = $ensp_index[prot] || $ensp_index_all[prot]
      [ensp, change]*":"
    end
    tsv = tsv.reorder "Mutated Isoform"
    tsv
  end

  dep :predict_ensp
  task :predict_all => :tsv do
    stream = TSV.get_stream step(:predict_ensp)
    s1, s2 = Misc.tee_stream stream
    FileUtils.mkdir_p files_dir
    smutations = CMD.cmd('cut -f 1', :in => s2, :pipe => true)
    smutations1, smutations = Misc.tee_stream smutations
    threads = []

    threads << Thread.new do 
      io = TSV.get_stream(Structure.job(:mi_interfaces, clean_name, :mutated_isoforms => smutations1).run(true))
      Open.write(file('interfaces'), io)
    end

    smutations2, smutations = Misc.tee_stream smutations
    threads << Thread.new do 
      io = TSV.get_stream(DbNSFP.job(:annotate, clean_name, :mutations => smutations2).run(true))
      Open.write(file('dbNSFP'), io)
    end

    databases = Structure::ANNOTATORS.keys - ['variants']

    databases.each do |database|
      smutations3, smutations = Misc.tee_stream smutations
      threads << Thread.new do 
        io = TSV.get_stream(Structure.job(:annotate_mi, clean_name + ': ' + database, :database => database, :mutated_isoforms => smutations3).run(true))
        Open.write(file(database), io)
      end
    end

    databases.each do |database|
      smutations3, smutations = Misc.tee_stream smutations
      threads << Thread.new do 
        io = TSV.get_stream(Structure.job(:annotate_mi_neighbours, clean_name + ': ' + database, :database => database, :mutated_isoforms => smutations3).run(true))
        Open.write(file(database + '_neighbours'), io)
      end
    end
    Misc.consume_stream smutations, true

    TmpFile.with_file() do |tmp_path|
      threads << Thread.new{Open.write(tmp_path, s1) }
      threads.each{|t| t.join }
      FileUtils.mv tmp_path, path
    end

    nil
  end

  export_asynchronous :predict, :predict_fix, :predict_all

end
