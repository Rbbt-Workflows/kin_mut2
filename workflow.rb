require 'rbbt/workflow'
require 'rbbt/util/misc/exceptions'

Workflow.require_workflow "Genomics"
module KinMut2
  extend Workflow

  helper :predict_script do
    Rbbt.bin["predict_mutations.pl"].find(:lib)
  end

  def self.organism
    Organism.default_code("Hsa")
  end


  input :mutations, :array, "Mutations"
  task :translate => :array do |mutations|
    raise ParameterException, "No mutations provided" if mutations.nil? or mutations.empty?

    translated = []
    organism = Organism.default_code("Hsa")
    translations = TSV.setup({}, :key_field => "Protein", :fields => ["UniProt/SwissProt Accession"], :type => :single, :namespace => organism)

    ensp = Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :fields => ["Ensembl Protein ID"], :persist => true)
    name = Organism.identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true)
    missing = []
    TSV.traverse mutations, :into => translated do |mutation|
      mutation = mutation.first if Array === mutation
      protein, change = mutation.split(/,|\s|:/)
      next unless change =~ /^[A-Z]\d+[A-Z]$/
      translation = case protein
                    when /ENSP/
                      translations[protein] = ensp[protein]
                    else
                      translations[protein] = name[protein] || protein
                    end
      if translation.nil?
        missing << protein
        next
      end
      [translation, change] * " "
    end

    Open.write(file(:translations), translations.to_s)
    set_info :missing_translations, missing

    translated
  end

  dep :translate
  task :predict => :tsv do 
    mutations = step(:translate).load
    iif mutations
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
    index = TSV.index(file(:translations), :target => "Protein")
    tsv = tsv.add_field "Fixed mutation" do |mutation,values|
      uni, change = mutation.split(" ")
      orig = index[uni]
      [orig, change]  * ":"
    end
    
    tsv = tsv.reorder("Fixed mutation", ["Prediction", "Score"])
    tsv.key_field = "Mutated Isoform"
    Open.write(file(:fixed), tsv.to_s)

    file_txt
  end

  dep :predict
  task :predict_fix => :tsv do
    TSV.get_stream step(:predict).file(:fixed)
  end

  export_asynchronous :predict, :predict_fix

end
