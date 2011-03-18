#!/usr/bin/perl
use strict;
use tgi_basic;
use tgi_bio;

my $data_path="/home/jmgonzalez/PROYECTOS/META_PATHOGENICITY_PREDS/WWW/data";
my $uniprot_path="$data_path/UNIPROT";
my @query=tgi_open($ARGV[0]);

my @numbered_list=tgi_open($ARGV[1]);
my %feature2index;
my %index2feature;
foreach my $lin (@numbered_list) {
	my ($n,$feat)=split(/\s+/,$lin);
	$feature2index{$feat}=$n;
	$index2feature{$n}=$feat;
	}

#Load the valid kinase accessions
my %kinase_accessions;
my @kinase_accessions=tgi_open("$data_path/KinaseAccessions.list");
foreach my $acc (@kinase_accessions) {
	$kinase_accessions{$acc}=1;
	}

#Get the accessions in the query
my %accs;
my @filtered_query;
foreach my $doublet (@query) {
	my ($acc,$wt,$mt,$pos,$mutationmention)=tgi_parse_doublet($doublet);
	$accs{$acc}=1;
	
	#Filter out non kinase proteins
	if ($kinase_accessions{$acc}!=1) {
		print STDERR "$acc is not a valid kinase identifier.... skipping\n";
		}
	push(@filtered_query,$doublet);
	}
my @accs=keys %accs;

#Kyte-Doolittle hydrophobicity
my %KDhydrophobicity;
load_hydrophobicity();

#Load kinase groups and sequences
my %uniprot_groups;
my %uniprot_sequences;
load_kinase_groups_and_sequences();

#Load Phosphoelm
my %phosphoelm;
tgi_load_phosphoelm();

#Load sumGOLOR
my %sumGOLOR;
load_sumGOLOR();

#Load UNIPROT/SWISSPROT annotations
my %uniprot_annotations;
tgi_load_swissprot_annotations(\@accs,\%uniprot_annotations);

#Load FireDB annotations
my %FireDB;
tgi_load_firedb();

#Load SIFT predictions (if exists reads from cache, otherwise queries the SIFT server)
my %SIFT_predictions;
my $sift_file="$data_path/SIFT.precomputed.txt";
tgi_load_SIFTwww_predictions(\@filtered_query,$sift_file,\%SIFT_predictions);

#Load PFAM domain information
my %PFAMdoms;
tgi_load_PFAMdomains(\@accs,\%PFAMdoms);

my %clusters;
my %acc2PFAM;
my $TDclusters="$data_path/TREEDETERMINANTS/Selected_Pfams_with_S3detOUTClusters";
tgi_TreeDeterminants_store_PFAM_clusters(\%clusters,\%acc2PFAM);
my %MSAs;
my $TDmsas="$data_path/TREEDETERMINANTS/Selected_Pfams_with_S3detOUTClusters_MulAlignments";
tgi_TreeDeterminants_store_MSAs(\%MSAs);

#Generate the vectors
my %vector;
foreach my $doublet (@filtered_query) {
	#Create a repository of accessions and doublets
	my ($acc,$wt,$mt,$pos,$mutationmention)=tgi_parse_doublet($doublet);
	my $acc_pos=$acc."_".$pos;
	
	#Reset the values of the features
	$vector{$doublet}{'TRAINING_LABEL'}=0;
	foreach my $feat (keys %feature2index) {
		$vector{$doublet}{$feat}=0;
		}
	$vector{$doublet}{'SIFTscore'}=0.0749; #Same value as undefined (SIFTwww) 
	$vector{$doublet}{'SIFTscore_binned'}=-1;
		
	#Add the AATYPE
	$vector{$doublet}{"aatype$wt"}=-1;
	$vector{$doublet}{"aatype$mt"}=1;
	
	#Add the HYDROPHOBICITY
	$vector{$doublet}{'biochem_diffkdhydrophobicity'}=$KDhydrophobicity{$mt}-$KDhydrophobicity{$wt};
	
	#Add the UNIPROT_GROUP
	$vector{$doublet}{"class_uniprotgroup_$uniprot_groups{$acc}"}=1;
	
	#Add PhosphoELM
	$vector{$doublet}{'phosphoelm'}=-1; #Not phosphorylable residues
	if ($phosphoelm{$acc_pos}==1) {
		$vector{$doublet}{'phosphoelm'}=1; #Phosphorylated Ser,Thr or Tyr
		}
	if (($phosphoelm{$acc_pos}!=1) && ($wt=~/[STY]/)) {
		$vector{$doublet}{'phosphoelm'}=0; #Not phosphorylated Ser,Thr or Tyr. 
		}
		
	#Add sumGOLOR
	$vector{$doublet}{'sumGOLOR'}=$sumGOLOR{$acc};
	
	#Add the Uniprot Annotations
	if ($uniprot_annotations{$acc}{'ACT_SITE'}{$pos}==1) {
		$vector{$doublet}{'swannot_act_site'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'CARBOHYD'}{$pos}==1) {
		$vector{$doublet}{'swannot_carbohyd'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'BINDING'}{$pos}==1) {
		$vector{$doublet}{'swannot_binding'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'CATALYTIC'}{$pos}==1) {
		$vector{$doublet}{'swannot_catalytic'}=1; 
		$vector{$doublet}{'swannot_any'}=1;
		}
	if ($uniprot_annotations{$acc}{'DISULFID'}{$pos}==1) {
		$vector{$doublet}{'swannot_disulfid'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'METAL'}{$pos}==1) {
		$vector{$doublet}{'swannot_metal'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'MOD_RES'}{$pos}==1) {
		$vector{$doublet}{'swannot_mod_res'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'MUTAGEN'}{$pos}==1) {
		$vector{$doublet}{'swannot_mutagen'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'NP_BIND'}{$pos}==1) {
		$vector{$doublet}{'swannot_np_bind'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'PTM'}{$pos}==1) {
		$vector{$doublet}{'swannot_ptm'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'SIGNAL'}{$pos}==1) {
		$vector{$doublet}{'swannot_signal'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}	
	if ($uniprot_annotations{$acc}{'SITE'}{$pos}==1) {
		$vector{$doublet}{'swannot_site'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	if ($uniprot_annotations{$acc}{'TRANSMEM'}{$pos}==1) {
		$vector{$doublet}{'swannot_transmem'}=1; 
		$vector{$doublet}{'swannot_any'}=1; 
		}
	
	#Add FireDB
	if ($FireDB{$acc_pos}==1) {
		$vector{$doublet}{'firedb'}=1; 
		}
		
	#Add SIFT score
	$vector{$doublet}{'SIFTscore'}=$SIFT_predictions{$doublet}{'SIFTscore'};
	$vector{$doublet}{'SIFTscore_binned'}=tgi_convert_SIFTscore($SIFT_predictions{$doublet}{'SIFTscore'}); 
	
	#Add PFAM domains
	foreach my $pfam_acc (sort keys %PFAMdoms) {
		if ($PFAMdoms{$pfam_acc}{$acc}{$pos}==1) {
			my $tag="pfam_".$pfam_acc;
			$vector{$doublet}{$tag}=1;
			}
		}
	}

#Add Tree Determinants separatedly
tgi_add_TreeDeterminants_to_training_vector_fscorebased(\@filtered_query,\%vector);



#Print the vector of features
foreach my $doublet (@filtered_query) {
	print "$vector{$doublet}{'TRAINING_LABEL'}";
	for (my $i=1;$i<=scalar keys %index2feature;$i++) {
		print "\t$i:$vector{$doublet}{$index2feature{$i}}";
		}
	print "\t#$doublet\n";
	}


#Exit properly
exit;

#####################################SUBROUTINES#####################################################

sub tgi_add_TreeDeterminants_to_training_vector_fscorebased {
	my $all_mutations=$_[0]; #Reference to the array that stores the doublets
	my $vector=$_[1]; #Reference to the hash that stores the features

	foreach my $doublet (@$all_mutations) {
		my ($acc,$wt,$mt,$pos,$mutationmention)=tgi_parse_doublet($doublet);
		$$vector{$doublet}{'TDs_fscore_wt'}=0;
		$$vector{$doublet}{'TDs_fscore_mt'}=0;
		$$vector{$doublet}{'TDs_fscore_diff'}=0;
		foreach my $PFAMdom (sort keys %{$acc2PFAM{$acc}}) {
			foreach my $acc_long_doublet (sort keys %{$acc2PFAM{$acc}{$PFAMdom}}) {
				my ($acc_again,$id,$organism,$start,$end,$length)=tgi_parse_acclong($acc_long_doublet);
				if (($pos>=$start) && ($pos<=$end)) {
					#Get the position of the mutation in the PFAM domain sequence (trimmed to the domain limits)
					my $pfam_pos=$pos-$start+1;
					#Get the position of the mutation in the MSA
					my $pfam_msaseq=$MSAs{$PFAMdom}{$acc_long_doublet}{'MSA'};
					my $pfam_msapos=tgi_get_MSA_equivalent_position($pfam_msaseq,$pfam_pos);

					#Calculate the conservation of the position for the rest of the cluster - wild type
					my $TPwt=0;
					my $TNwt=0;
					my $FPwt=0;
					my $FNwt=0;
					foreach my $acc_long_msa (sort keys %{$MSAs{$PFAMdom}}) {
						next if ($acc_long_msa eq $acc_long_doublet);
						my $msa_seq=$MSAs{$PFAMdom}{$acc_long_msa}{'MSA'};
						my @msa_seq=split("",$msa_seq);
						my $residue_msa=$msa_seq[$pfam_msapos];
						if ($clusters{$PFAMdom}{$acc_long_msa}{$acc_long_doublet}==1) { #Same cluster/subfamily
							if ($wt eq $residue_msa) {
								$TPwt++; #True positive: Same cluster, same residue
								}
							if ($wt ne $residue_msa){
								$FNwt++; #False negative: Same cluster, different residue
								}
							}
						if ($clusters{$PFAMdom}{$acc_long_msa}{$acc_long_doublet}!=1) {
							if ($wt eq $residue_msa) {
								$FPwt++; #False positive: Different cluster, same residue
								}
							if ($wt ne $residue_msa){
								$TNwt++; #True negative: Different cluster, different residue
								}
							}
						}
					my ($accuracy_wt,$recall_wt,$precision_wt,$mcc_wt,$fscore_wt,$n_wt)=tgi_evaluate_performance($TPwt,$TNwt,$FPwt,$FNwt);

					#Calculate the conservation of the position for the rest of the cluster - mutant
					my $TPmt=0;
					my $TNmt=0;
					my $FPmt=0;
					my $FNmt=0;
					foreach my $acc_long_msa (sort keys %{$MSAs{$PFAMdom}}) {
						next if ($acc_long_msa eq $acc_long_doublet);
						my $msa_seq=$MSAs{$PFAMdom}{$acc_long_msa}{'MSA'};
						my @msa_seq=split("",$msa_seq);
						my $residue_msa=$msa_seq[$pfam_msapos];
						if ($clusters{$PFAMdom}{$acc_long_msa}{$acc_long_doublet}==1) { #Same cluster/subfamily
							if ($mt eq $residue_msa) {
								$TPmt++; #True positive: Same cluster, same residue
								}
							if ($mt ne $residue_msa){
								$FNmt++; #False negative: Same cluster, different residue
								}
							}
						if ($clusters{$PFAMdom}{$acc_long_msa}{$acc_long_doublet}!=1) {
							if ($mt eq $residue_msa) {
								$FPmt++; #False positive: Different cluster, same residue
								}
							if ($mt ne $residue_msa){
								$TNmt++; #True negative: Different cluster, different residue
								}
							}
						}
					my ($accuracy_mt,$recall_mt,$precision_mt,$mcc_mt,$fscore_mt,$n_mt)=tgi_evaluate_performance($TPmt,$TNmt,$FPmt,$FNmt);
					my $diff_fscore=$fscore_mt-$fscore_wt;

					#Add the results to the training vector
					$$vector{$doublet}{'TDs_fscore_wt'}=$fscore_wt;
					$$vector{$doublet}{'TDs_fscore_mt'}=$fscore_mt;
					$$vector{$doublet}{'TDs_fscore_diff'}=$diff_fscore;
					}
				}
			}
		}
	return 0;
	}


sub tgi_TreeDeterminants_store_PFAM_clusters {
	my $clusters=$_[0]; #Reference to the hash that will store the clusters
	my $acc2PFAM=$_[1]; #Reference to the hash that will store the acc<=>PFAM equivalence
	my @PFAMdoms_withcluster=`cd $TDclusters && ls PF*`;
	foreach my $PFAMdom_withcluster_file (@PFAMdoms_withcluster) {
		chomp($PFAMdom_withcluster_file);
		#Get the PFAM identifier from the filename
		my $PFAMdom=tgi_getPFAMfromfilename($PFAMdom_withcluster_file);
		#Store the file content
		my @PFAMdom_withcluster_file=tgi_open("$TDclusters/$PFAMdom_withcluster_file");
		my @cluster_lines=grep (/^CL:\s+/,@PFAMdom_withcluster_file);
		for (my $i;$i<scalar @cluster_lines-1;$i++) {
			my ($null,$acc_long1,$cluster1)=split(/\s+/,$cluster_lines[$i]);
			$$clusters{$PFAMdom}{$acc_long1}{$acc_long1}=1;
			my $acc1=substr($acc_long1,0,6);
			$$acc2PFAM{$acc1}{$PFAMdom}{$acc_long1}=1;
			#print STDERR "Storing $acc1 - $PFAMdom - $acc_long1\n";
			for (my $j=$i+1;$j<scalar @cluster_lines;$j++) {
				my ($null,$acc_long2,$cluster2)=split(/\s+/,$cluster_lines[$j]);
				$$clusters{$PFAMdom}{$acc_long2}{$acc_long2}=1;
				my $acc2=substr($acc_long2,0,6);
				$$acc2PFAM{$acc2}{$PFAMdom}{$acc_long2}=1;
				#print STDERR "Storing $acc2 - $PFAMdom - $acc_long2\n";
				if ($cluster1==$cluster2) {
					$$clusters{$PFAMdom}{$acc_long1}{$acc_long2}=1;
					$$clusters{$PFAMdom}{$acc_long2}{$acc_long1}=1;
					}
				}
			}
		}
	return 0;
	}

sub tgi_TreeDeterminants_store_MSAs {
	my $MSAs=$_[0]; #Reference to the hash that will store the MSAs
	my @PFAMdoms_msas=`cd $TDmsas && ls PF*`;
	foreach my $PFAMmsa_file (@PFAMdoms_msas) {
		chomp($PFAMmsa_file);
		#Get the PFAM identifier from the filename
		my $PFAMdom=tgi_getPFAMfromfilename($PFAMmsa_file);
		#Store the file content
		my @PFAMmsa_file=tgi_open("$TDmsas/$PFAMmsa_file");
		foreach my $lin (@PFAMmsa_file) {
			my ($acc_long,$msa_seq)=split(/\s+/,$lin);
			my $seq=$msa_seq;
			$seq=~s/-//g;
			$$MSAs{$PFAMdom}{$acc_long}{'SEQ'}=$seq;
			$$MSAs{$PFAMdom}{$acc_long}{'MSA'}=$msa_seq;
			}
		}
	return 0;
	}

sub tgi_parse_acclong {
	my $acc_long=$_[0];
	my $acc=substr($acc_long,0,6);
	my ($null,$rest)=split("=",$acc_long);
	my ($id,$pos)=split("/",$rest);
	my ($null,$organism)=split("_",$id);
	my ($start,$end)=split("-",$pos);
	my $length=$end-$start+1;
	return ($acc,$id,$organism,$start,$end,$length)
	}

sub tgi_getPFAMfromfilename {
	my $filename=$_[0];
	$filename=~/(PF\d\d\d\d\d)/;
	my $pfam=$1;
	return $pfam;
	}

sub tgi_get_MSA_equivalent_position($$) {
	my $pfam_msaseq=$_[0];
	my $target_pfam_pos=$_[1];
	my $msapos=0;
	my @pfam_msaseq=("X",split("",$pfam_msaseq));#So the position in the seq and the position in the array is teh same
	for (my $i=1;$i<scalar @pfam_msaseq;$i++) {
		$msapos++ if ($pfam_msaseq[$i] ne "-");
		return $i if ($msapos == $target_pfam_pos);
		}
	return 0;
	}

sub tgi_load_SIFTwww_predictions {
	use URI::Escape;
	my $all_doublets=$_[0]; #Reference to the array with all the mutations
	my $sift_file=$_[1]; #Path to the file with the pre-computed SIFT scores
	my $sift=$_[2]; #Reference to the hash that will store the information

	#If the SIFT_FILE does not exist, create one
	unless (-e $sift_file) {
		system("touch $sift_file");
		}

	#Load the precomputed SIFT scores
	my @precomputed=tgi_open($sift_file);
	my %precomputed;
	foreach my $lin (@precomputed) {
		my ($doublet,$np,$SIFTscore,$n_seqs)=split("\t",$lin);
		$$sift{$doublet}{'SIFTscore'}=$SIFTscore;
		$$sift{$doublet}{'SIFTnseqs'}=$n_seqs;
		$precomputed{$doublet}=1;
		}

	my %buffer;
	my %acc2np;
	my %np2acc; 

	#buffer the query
	foreach my $doublet (@$all_doublets) {
		if (defined $$sift{$doublet}{'SIFTscore'}) {
			next;
			}
		#Store the "undefined" values, just in case
		$$sift{$doublet}{'SIFTscore'}=0.0749;
		$$sift{$doublet}{'SIFTnseqs'}=0;
		my ($acc,$wt,$mt,$pos,$mutationmention)=tgi_parse_doublet($doublet);

		#Get the NP_ from the accession
		if (!defined $acc2np{$acc}) {
			my @sfetch=grep (/DR\s*RefSeq;\s*NP_/,`sfetch -Dsw $acc`);
			$sfetch[0]=~/(NP_\d+)/;
			$acc2np{$acc}=$1;
			$np2acc{$1}=$acc;
			}
		my $np=$acc2np{$acc};
		$buffer{$acc}{$mutationmention}=1;
		}

	#Run the query
	my $ct=0;
	my $url_buffer;
	my $query;
	my %checked_accs;
	my @buffer=sort keys %buffer;
	for (my $i=0;$i<scalar @buffer;$i++) {
		my $acc=$buffer[$i];
		next if ($checked_accs{$acc}==1);
		next if (!defined $acc2np{$acc});
		$ct++;
		$checked_accs{$acc}=1;
		my $np=$acc2np{$acc};
		my $muts=join(",",sort keys %{$buffer{$acc}});
		$url_buffer.="$np,$muts\n";
		$query.="$np\t$acc\t$muts\n";
		if (($ct==50) || ($i==scalar @buffer-1)) {
			
			#Run the query
			chomp($url_buffer);
			my $buffer_encoded=uri_escape($url_buffer);
			my @sift_html=`wget http://sift.jcvi.org/sift-bin/SIFT_pid_subst_all_submit.pl -O - --post-data=\"GI=$buffer_encoded&sequences_to_select=BEST&seq_identity_filter=90\"`;

			#Parse the query
			my @np_lines=grep (/^<td>NP_/,@sift_html);
			foreach my $results (@np_lines) {
				chomp($results);
				$results=~s/^<td>//gi;
				$results=~s/<\/tr>$//gi;
				$results=~s/<\/font>//gi;
				$results=~s/<\/td><td>/\t/gi;
				$results=~s/<font color=red>//gi;
				$results=~s/<\/td>//gi;
				my ($np,$mutationmention,$dbsnp,$prediction,$siftscore,$median,$n_seqs)=split("\t",$results);
				my $acc=$np2acc{$np};
				my $doublet=$acc."_".$mutationmention;

				$$sift{$doublet}{'SIFTscore'}=$siftscore;
				$$sift{$doublet}{'SIFTnseqs'}=$n_seqs;

				#Exceptions: Less than 10 sequences, very neutral score.
				if ($$sift{$doublet}{'SIFTnseqs'}<=10) {
					$$sift{$doublet}{'SIFTscore'}=0.0749;
					}

				#Store the result to be precomputed in the future
				$buffer{$acc}{$mutationmention}=0;
				next if ($precomputed{$doublet}==1);
				open (SIFT,">>$sift_file");
				print SIFT "$doublet\t$np\t$$sift{$doublet}{'SIFTscore'}\t$$sift{$doublet}{'SIFTnseqs'}\n";
				close (SIFT);
				}

			#Reset values
			$ct=0;
			$url_buffer="";
			$query="";
			}
		}
	
	#Force the results that were requiered (buffered) but not retrieved somehow in the analysis. 
	foreach my $acc (sort keys %buffer) {
		foreach my $mutationmention (sort keys %{$buffer{$acc}}) {
			if ($buffer{$acc}{$mutationmention}==1) { #Buffered but not retrieved in tha analysis
				my $doublet=$acc."_".$mutationmention;
				my $np=$acc2np{$acc};
				#print STDERR "Forcing $doublet ($np)\n";
				$$sift{$doublet}{'SIFTscore'}=10;
				$$sift{$doublet}{'SIFTnseqs'}=0;
				open (SIFT,">>$sift_file");
				print SIFT "$doublet\t$np\t$$sift{$doublet}{'SIFTscore'}\t$$sift{$doublet}{'SIFTnseqs'}\n";
				close (SIFT);
				}
			}
		}
	
	return 0;
	}

sub tgi_load_firedb {
	my @firedb=tgi_open("$data_path/FireDB.txt");
	foreach my $lin (@firedb) {
		my @fields=split("\t",$lin);
		my $acc=$fields[0];
		for (my $i=1;$i<scalar @fields;$i++) {
			my $acc_pos=$acc."_".$fields[$i];
			$FireDB{$acc_pos}=1;
			}
		}
	}

sub load_sumGOLOR {
	my @sumGOLOR=tgi_open("$data_path/sumGOLOR.txt");
	foreach my $lin (@sumGOLOR) {
		my ($acc,$type,$sumGOLOR)=split(/\s+/,$lin);
		$sumGOLOR{$acc}=$sumGOLOR;
		}
	return 0;
	}

sub tgi_load_phosphoelm {
	my @phosphoelm=tgi_open("$data_path/phosphoelm.txt");
	foreach my $lin (@phosphoelm) {
		next if ($lin=~/^#/);
		my ($acc,$pos,$aa,$validation)=split("\t",$lin);
		next if ($validation ne "X");
		my $tag=$acc."_".$pos;
		$phosphoelm{$tag}=1;
		}
	return 0;
	}

sub load_kinase_groups_and_sequences {
	my @basic_info=tgi_open("$data_path/KinaseAccessions_Group_Seqs.txt");
	foreach my $lin (@basic_info) {
		my ($acc,$group,$seq)=split("\t",$lin);
		$uniprot_groups{$acc}=$group;
		$uniprot_sequences{$acc}{$seq};
		}
	return 0;
	}

sub tgi_evaluate_performance {
	my $TP=$_[0];
	my $TN=$_[1];
	my $FP=$_[2];
	my $FN=$_[3];
	my $n=($TP+$TN+$FN+$FP);
	my $accuracy=100*($TP+$TN)/($TP+$TN+$FN+$FP);
	my $precision=0;
	if (($TP+$FP)!=0) {
		$precision=100*($TP)/($TP+$FP);
		}
	my $recall=0;
	if (($TP+$FN)!=0) {
		$recall=100*($TP)/($TP+$FN);
		}
	my $mcc=0;
	if (($TP+$FN)*($TP+$FN)*($TN+$FN)*($TN+$FP)!=0) {
		$mcc=(($TP*$TN)-($FP*$FN))/sqrt(($TP+$FN)*($TP+$FN)*($TN+$FN)*($TN+$FP));
		}
	my $fscore=0;
	if ($precision*$recall!=0) {
		$fscore=2/((1/$recall)+(1/$precision));
		}
	return ($accuracy,$recall,$precision,$mcc,$fscore,$n);
	}

sub load_hydrophobicity {
	$KDhydrophobicity{'I'}=4.5;
	$KDhydrophobicity{'V'}=4.2;
	$KDhydrophobicity{'L'}=3.8;
	$KDhydrophobicity{'F'}=2.8;
	$KDhydrophobicity{'C'}=2.5;
	$KDhydrophobicity{'M'}=1.9;
	$KDhydrophobicity{'A'}=1.8;
	$KDhydrophobicity{'G'}=-0.4;
	$KDhydrophobicity{'T'}=-0.7;
	$KDhydrophobicity{'S'}=-0.8;
	$KDhydrophobicity{'W'}=-0.9;
	$KDhydrophobicity{'Y'}=-1.3;
	$KDhydrophobicity{'P'}=-1.6;
	$KDhydrophobicity{'H'}=-3.2;
	$KDhydrophobicity{'E'}=-3.5;
	$KDhydrophobicity{'Q'}=-3.5;
	$KDhydrophobicity{'N'}=-3.5;
	$KDhydrophobicity{'D'}=-3.5;
	$KDhydrophobicity{'K'}=-3.9;
	$KDhydrophobicity{'R'}=-4.5;
	}

sub tgi_open {
	my $file=$_[0];
	open (TGI_OPEN,"<$file") or die("ERROR: Can't open: $file\n") ;
	chomp(my @file=<TGI_OPEN>);
	close(TGI_OPEN);
	return @file;
	}

sub tgi_parse_doublet {
	my $doublet=$_[0];
	my ($acc,$mutationmention)=split("_",$doublet);
	my @split=split("",$mutationmention);
	my $wt=$split[0];
	my $mt=$split[scalar @split-1];
	my $pos;
	for (my $i=1;$i<scalar @split-1;$i++) {
		$pos.=$split[$i];
		}
	return ($acc,$wt,$mt,$pos,$mutationmention);
	}

sub tgi_load_swissprot_annotations {
	my $kinase_accs=$_[0]; #Reference to the array with the accessions
	my $uniprot_annotations=$_[1]; #Reference to the array that will store the annotations

	foreach my $acc (@$kinase_accs) {
		my $uniprot_file="$uniprot_path/uniprot.$acc.txt";
		unless (-e $uniprot_file) {
			my @sfetch=`sfetch -Dsw $acc > $uniprot_file`;
			}
		my @uniprot_file=tgi_open($uniprot_file);
		foreach my $ft_lin (grep (/^FT/,@uniprot_file)) {
			my ($ft,$type,$start,$end)=split(/\s+/,$ft_lin);
			if ($type eq 'SITE') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'SITE'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'BINDING') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'BINDING'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'ACT_SITE') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'ACT_SITE'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'LIPID') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'LIPID'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'METAL') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'METAL'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'CARBOHYD') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'CARBOHYD'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'DNA_BIND') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'DNA_BIND'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'NP_BIND') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'NP_BIND'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'CA_BIND') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'CA_BIND'}{$i}=1;
					$$uniprot_annotations{$acc}{'CATALYTIC'}{$i}=1;
					}
				next;
				}
			if ($type eq 'DISULFID') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'DISULFID'}{$i}=1;
					}
				next;
				}
			if ($type eq 'MOD_RES') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'MOD_RES'}{$i}=1;
					$$uniprot_annotations{$acc}{'PTM'}{$i}=1;
					}
				next;
				}
			if ($type eq 'PROPEP') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'PROPEP'}{$i}=1;
					$$uniprot_annotations{$acc}{'PTM'}{$i}=1;
					}
				next;
				}
			if ($type eq 'SIGNAL') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'SIGNAL'}{$i}=1;
					$$uniprot_annotations{$acc}{'PTM'}{$i}=1;
					}
				next;
				}
			if ($type eq 'MUTAGEN') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'MUTAGEN'}{$i}=1;
					}
				next;
				}
			if ($type eq 'TRANSMEM') {
				for (my $i=$start;$i<=$end;$i++) {
					$$uniprot_annotations{$acc}{'TRANSMEM'}{$i}=1;
					}
				next;
				}
			}
		}
	return 0;
	}

sub tgi_convert_SIFTscore {
	my $raw_SIFTscore=$_[0];
	if ($raw_SIFTscore==-1) {
		return 0; #Undefined
		}
	if ($raw_SIFTscore<=0.05) {
		return 1; #Pathogenic
		}
	if ($raw_SIFTscore>0.05) {
		return -1; #Neutral
		}
	return 0;
	}

sub tgi_load_PFAMdomains {
	my $kinase_accs=$_[0]; #Reference to the array with the accessions
	my $PFAMdoms=$_[1]; #Reference to the hash that will store the information

	#Store the kinases in a hash
	my %kinase_accs;
	foreach my $acc (@$kinase_accs) {
		$kinase_accs{$acc}=1;
		}

	#Parse the domains file
	my $acc;
	open (PFAMFILE,"<$data_path/swisspfam.kinasefiltered");
	while (my $lin=<PFAMFILE>) {
		chomp($lin);

		#Empty lines
		next if ($lin=~/^\s*$/); 

		#Protein lines
		if ($lin=~/^>/) {
			my @fields=split(/\s+/,$lin);
			$acc=$fields[2];
			$acc=~s/\.\d+//;
			}

		next if ($kinase_accs{$acc}!=1); #Skip lines that do not correspond to kinases

		#Domain definition lines	
		if ($lin=~/\)\s*(PF.*)/) { #Pfam A only
			my $info=$1;
			my @fields=split(/\s+/,$info);
			my $pfam_acc=$fields[0];
			$pfam_acc=~s/\.\d+//;


			#Record the domain positions
			for (my $i=1;$i<scalar @fields;$i++) {
				next if ($fields[$i]!~/\d+-\d+/);
				my ($start,$end)=split("-",$fields[$i]);
				for (my $pos=$start;$pos<=$end;$pos++) {
					$$PFAMdoms{$pfam_acc}{$acc}{$pos}=1;
					$$PFAMdoms{'any'}{$acc}{$pos}=1;
					}
				}
			}
		}
	close(PFAMFILE);
	return 0;
	}