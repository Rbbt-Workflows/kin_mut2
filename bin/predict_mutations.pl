#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

#Define the default settings of the script here
my $version_tag="1";
my $version=0;
my $help=0;
my $step_ct=0;
my $prediction_input="";
my $output_path="";

my $base_dir= dirname(dirname(__FILE__));
my $uniprot_repository="$base_dir/data/UNIPROT";
my $kinase_basicinfo="$base_dir/data/Uniprot.kinase.basicinfo.txt";
my $kinase_domains_lor="$base_dir/data/Uniprot.kinase.domains_lor.txt";
my $kinase_groups_lor="$base_dir/data/Uniprot.kinase.groups_lor.txt";
my $kinase_sumgolor="$base_dir/data/Gene_Ontology.kinase.sumGOLOR.txt";

#Built-in configuration
my $weka_jar="$base_dir/share/software/opt/jar/weka.jar";
my $weka_random_forest="$base_dir/share/model/weka.randomforest.model";
my $sift_na="0.0749";


#Collect the input from the command line, force exit if minimun variables not specified, print help if needed
GetOptions (
		'input=s' => \$prediction_input,
		'output=s' => \$output_path,
		'help' => \$help,
		'version'=> \$version
		);	

if ($prediction_input eq "") {
	print STDERR "\nERROR: Input file (mutations to predict) must be specified\n";
	tgi_print_help();
	}
unless (-e $prediction_input) {
	print STDERR "\nERROR: Input file (mutations to predict) not found: $prediction_input\n";
	tgi_print_help();
	}		
if ($output_path eq "") {
	print STDERR "\nERROR: Destination folder must be specified\n";
	tgi_print_help();
	}
tgi_create_dir($output_path);

tgi_print_help() if ($help);
tgi_print_version() if ($version) ;

my $weka_vector_file="$output_path/vectors.weka.arff";
my $weka_vectors_out="$output_path/vectors.weka.out";
my $weka_predictions="$output_path/vectors.weka.predictions";

print STDERR "#PERL SCRIPT:\t$0\n";
print STDERR "#\n";

#Load the biochemical properties of the amino acids
my %volumes;
my %cbetabranching;
my %KDhydrophobicity;
my %formalcharge;
tgi_biochemical_properties();

my %kinase_features;
my %kinase_position_features;

#Parse the basic info file so that we get the equivalences and other information
$step_ct++;
print STDERR "#STEP $step_ct: Parsing basic information file\n";
my @basic=tgi_open($kinase_basicinfo);
for (my $i=1;$i<scalar @basic;$i++) {
	my ($gene,$id,$acc,$ensp,$group)=split(/\t/,$basic[$i]);
	$kinase_features{$acc}{'GROUP'}=$group;
	$kinase_features{$acc}{'ENSP'}=$ensp;
	$kinase_features{$acc}{'UNIPROT_ID'}=$id;
	}
undef(@basic);


##################################### KINASE-LEVEL FEATURES ###################################
#Get the sumGOLOR for each of the kinases
$step_ct++;
print STDERR "#STEP $step_ct: Parsing sumGOLOR precalculations\n";
my @sumGOLOR=tgi_open($kinase_sumgolor);
for (my $i=1;$i<scalar @sumGOLOR;$i++) {
	my ($acc,$sumgolor,$sumgolor_mf)=split(/\t/,$sumGOLOR[$i]);
	$kinase_features{$acc}{'SUMGOLOR'}=$sumgolor;
	$kinase_features{$acc}{'SUMGOLOR_MF'}=$sumgolor_mf;
	}
undef(@sumGOLOR);

#Get the log-odds ration at the GROUP level
$step_ct++;
print STDERR "#STEP $step_ct: Parsing group log-odds ratio precalculations\n";
my %grouplor;
my @grouplor=tgi_open($kinase_groups_lor);
for (my $i=1;$i<scalar @grouplor;$i++) {
	my @fields=split(/\t/,$grouplor[$i]);
	my $group=$fields[0];
	my $grouplor=$fields[5];
	$grouplor{$group}=$grouplor;
	}
undef(@grouplor);

##################################### DOMAIN-LEVEL FEATURES ###################################

#Get the log-odds ration at the DOMAIN level
$step_ct++;
print STDERR "#STEP $step_ct: Parsing domain log-odds ratio precalculations\n";
my %domlor;
my @domlor=tgi_open($kinase_domains_lor);
for (my $i=1;$i<scalar @domlor;$i++) {
	my @fields=split(/\t/,$domlor[$i]);
	my $domain=$fields[0];
	my $domlor=$fields[5];
	$domlor{$domain}=$domlor;
	}
undef(@domlor);


##################################### POSITION-LEVEL FEATURES ###################################

#Get the positional Uniprot, FireDB, Phosphoelm and domain features 
$step_ct++;
print STDERR "#STEP $step_ct: Obtaining the positional Uniprot, FireDB, Phosphoelm and domain features\n";
foreach my $acc (keys %kinase_features) {
	my @features=tgi_open("$uniprot_repository/$acc.features.txt");
	foreach my $lin (@features) {
		my ($pos,$feature,$desc)=split(/\t/,$lin);
		$kinase_position_features{$acc}{$pos}{$feature}=$desc;
		}
	}

#Get the positional SIFT scores
$step_ct++;
print STDERR "#STEP $step_ct: Obtaining the pathogenicity from SIFT\n";
my %SIFTscores;
foreach my $acc (keys %kinase_features) {
	my @sift=tgi_open("/$acc.pathogenicity.txt");
	foreach my $lin (@sift) {
		my ($acc,$mention,$sift)=split(/\t/,$lin); 
		$SIFTscores{$acc}{$mention}=$sift;
		}
	}

##################################### CREATE THE VECTORS ###################################
$step_ct++;
print STDERR "#STEP $step_ct: Generating the vectors for the complete dataset\n";
my %mutation_mentions;

my @prediction_input=tgi_open($prediction_input);
foreach my $lin (@prediction_input) {
	next if ($lin=~/^\s*$/);
	my ($acc,$mention)=split(/\s/,$lin);
	unless (defined $kinase_features{$acc}{'GROUP'}) {
		print STDERR "#\tWARNING: Skipping $lin - $acc is not a supported kinase accession\n";
		next;
		}
	$mutation_mentions{$acc}{$mention}{'group'}=$kinase_features{$acc}{'GROUP'};
	$mutation_mentions{$acc}{$mention}{'group_lor'}=$grouplor{$kinase_features{$acc}{'GROUP'}};
	$mutation_mentions{$acc}{$mention}{'type'}="?";
	}

my @uniprot_annotations=qw(ACT_SITE BINDING CARBOHYD DISULFID METAL MOD_RES NP_BIND REPEAT SIGNAL SITE TRANSMEM ZN_FING); 
my @features_to_print=qw(group group_lor sumGOlor domain domain_lor aa_wt aa_mt volume cbetabranching hydrophobicity formal_charge firedb phosphoelm any_uniprot act_site binding carbohyd disulfid metal mod_res np_bind repeat signal site transmem zn_fing sift type) ; ########################### DEBUG!!!!! 

foreach my $acc (sort keys %mutation_mentions) {
	foreach my $mention (keys %{$mutation_mentions{$acc}}) {
		
		my ($wt,$pos,$mt)=tgi_parse_mention($mention);
		$mutation_mentions{$acc}{$mention}{'aa_wt'}=$wt; # Wild type AA
		$mutation_mentions{$acc}{$mention}{'aa_mt'}=$mt; # Mutant AA
		
		$mutation_mentions{$acc}{$mention}{'volume'}=$volumes{$mt}-$volumes{$wt}; #Volume
		$mutation_mentions{$acc}{$mention}{'cbetabranching'}=$cbetabranching{$mt}-$cbetabranching{$wt}; #Cbeta Branching
		$mutation_mentions{$acc}{$mention}{'hydrophobicity'}=$KDhydrophobicity{$mt}-$KDhydrophobicity{$wt}; #Hydrophobicity
		$mutation_mentions{$acc}{$mention}{'formal_charge'}=$formalcharge{$mt}-$formalcharge{$wt}; #Formal charge
		
		$mutation_mentions{$acc}{$mention}{'sumGOlor'}=$kinase_features{$acc}{'SUMGOLOR'}; # sumGOlor
		
		$mutation_mentions{$acc}{$mention}{'domain'}='X';
		$mutation_mentions{$acc}{$mention}{'domain_lor'}=0;
		if (defined $kinase_position_features{$acc}{$pos}{'DOMAIN'}) {
			$mutation_mentions{$acc}{$mention}{'domain'}=$kinase_position_features{$acc}{$pos}{'DOMAIN'}; # Domain
			$mutation_mentions{$acc}{$mention}{'domain_lor'}=$domlor{$kinase_position_features{$acc}{$pos}{'DOMAIN'}}; #Domain lor
			}
		
		$mutation_mentions{$acc}{$mention}{'firedb'}=(defined ($kinase_position_features{$acc}{$pos}{'FIREDB'}) ? 1 : 0); #FireDB
		$mutation_mentions{$acc}{$mention}{'phosphoelm'}=(defined ($kinase_position_features{$acc}{$pos}{'PHOSPHOELM'}) ? 1 : 0); #PhosphoELM
		
		$mutation_mentions{$acc}{$mention}{'any_uniprot'}=0;
		foreach my $uniprot_annot (@uniprot_annotations) {
			my $uniprof_annot_lc=lc($uniprot_annot);
			$mutation_mentions{$acc}{$mention}{$uniprof_annot_lc}=(defined ($kinase_position_features{$acc}{$pos}{$uniprot_annot}) ? 1 : 0); #Uniprot Annotations
			$mutation_mentions{$acc}{$mention}{'any_uniprot'}=1 if (defined $kinase_position_features{$acc}{$pos}{$uniprot_annot});
			}
		
		$mutation_mentions{$acc}{$mention}{'sift'}=($SIFTscores{$acc}{$mention}!~/^\s*$/ ? $SIFTscores{$acc}{$mention} : $sift_na); #Sift
		}
	}

##################################### PRINT THE VECTORS IN WEKA FORMAT ###################################
$step_ct++;
print STDERR "#STEP $step_ct: Printing the vectors in WEKA format: $weka_vector_file\n";

my $groups=join(",",sort keys %grouplor);
my $aas=join(",", sort keys %volumes);

#\@attribute group {$groups}
open (WEKA, ">$weka_vector_file");
print WEKA <<EOF;
\@relation kinase_mutations_full_set
\@attribute group_lor real
\@attribute sumGOlor real
\@attribute domain_lor real
\@attribute aa_wt {$aas}
\@attribute aa_mt {$aas}
\@attribute volume real
\@attribute cbetabranching real
\@attribute hydrophobicity real
\@attribute formal_charge real
\@attribute firedb {TRUE,FALSE}
\@attribute phosphoelm {TRUE,FALSE}
\@attribute any_uniprot {TRUE,FALSE}
\@attribute act_site {TRUE,FALSE}
\@attribute binding {TRUE,FALSE}
\@attribute carbohyd {TRUE,FALSE}
\@attribute disulfid {TRUE,FALSE}
\@attribute metal {TRUE,FALSE}
\@attribute mod_res {TRUE,FALSE}
\@attribute np_bind {TRUE,FALSE}
\@attribute repeat {TRUE,FALSE}
\@attribute signal {TRUE,FALSE}
\@attribute site {TRUE,FALSE}
\@attribute transmem {TRUE,FALSE}
\@attribute zn_fing {TRUE,FALSE}
\@attribute sift real
\@attribute type {disease,neutral}
EOF
print WEKA "\n\n\@data\n";
my %instances;
my $instance_ct=0;
foreach my $acc (sort keys %mutation_mentions) {
	foreach my $mention (keys %{$mutation_mentions{$acc}}) {
		$instance_ct++;
		$instances{$instance_ct}="$acc $mention";
		foreach my $feat (@features_to_print) {
			my $value=$mutation_mentions{$acc}{$mention}{$feat};
			
			if (($feat ne "volume") && ($feat ne "formal_charge") && ($feat ne "cbetabranching") && ($feat ne "hydrophobicity") && ($feat!~/lor$/) && ($feat ne "sift")) {
				$value="TRUE" if ($value=~/^1$/);
				$value="FALSE" if ($value=~/^0$/);
				}
				
			next if ($feat eq "domain");
			next if ($feat eq "type");
			next if ($feat eq "group");
			
			$value="?" if ($value=~/^\s*$/);
 			print WEKA "$value,";
			}
		print WEKA "$mutation_mentions{$acc}{$mention}{'type'}\n";
		}
	}
close (WEKA);

################################################# PREDICT WITH WEKA #################################################################
$step_ct++;
print STDERR "#STEP $step_ct: Predicitons with trained random forest: $weka_random_forest\n";
print("java -Xmx1000M weka.classifiers.trees.RandomForest -l $weka_random_forest -T $weka_vector_file -p 0 > $weka_vectors_out");
system("java -classpath '$weka_jar' -Xmx1000M weka.classifiers.trees.RandomForest -l $weka_random_forest -T $weka_vector_file -p 0 > $weka_vectors_out");

open (PRED,">$weka_predictions");
my @weka_vectors_out=tgi_open($weka_vectors_out);
for (my $i=5;$i<scalar @weka_vectors_out;$i++) {
	next if ($weka_vectors_out[$i]=~/^\s*$/);
	$weka_vectors_out[$i]=~s/^\s+//;
	my @fields=split(/\s+/,$weka_vectors_out[$i]);
	my $instance=$fields[0];
	my $class=$fields[2];
	$class=~s/^\d://;
	my $score=$fields[3];
	print PRED "$instances{$instance}\t$class\t$score\n";
	}
close (PRED);

#Exit properly 
print STDERR "#\n#Bye!\n";
exit;



############################################### SUBROUTINES ####################################################
sub tgi_print_help 	{
	print "USAGE:\nperl $0 -input <INPUT_FILE> -output <DESTINATION_PATH>\n";
	print "OPTIONS:\n";
	print "\tinput\tinput file with the mutations to predict\n";
	print "\toutput\tdestination directory\n";
	print "\n";
	print "\tversion\tprints the version of the script\n";
	print "\thelp\tprints this help\n";
	exit;
	}

sub tgi_print_version {
	print "$0 Version: $version_tag\nJose MG Izarzugaza (txema\@cbs.dtu.dk)\n";
	exit;
	}

sub tgi_parse_mention {
	#H52R
	my $mention=$_[0];
	my $wt=uc(substr($mention,0,1));
	my $mt=uc(substr($mention,-1,1));
	my $pos=$mention;
	$pos=~s/[A-Z|a-z]//gi;

	return ($wt,$pos,$mt);
	}
	
sub tgi_open {
	my $file=$_[0];
	open (TGI_OPEN,"<$file") ;
	chomp(my @file=<TGI_OPEN>);
	close(TGI_OPEN);
	return @file;
	}

sub tgi_create_dir{
	system("mkdir -p $_[0]") unless (-e $_[0]) ;
	}

sub tgi_biochemical_properties {
	# USAGE:
	# my %volumes;
	# my %cbetabranching; #C-beta branching (number of non-hydrogen substituents attached to C-beta carbon)
	# my %KDhydrophobicity; #Kyte-Doolittle hydrophobicity
	# my %formalcharge;
	# tgi_biochemical_properties();

	$volumes{'A'}=88.6;
	$volumes{'R'}=173.4;
	$volumes{'D'}=111.1;
	$volumes{'N'}=114.1;
	$volumes{'C'}=108.5;
	$volumes{'E'}=138.4;
	$volumes{'Q'}=143.8;
	$volumes{'G'}=60.1;
	$volumes{'H'}=153.2;
	$volumes{'I'}=166.7;
	$volumes{'L'}=166.7;
	$volumes{'K'}=168.6;
	$volumes{'M'}=162.9;
	$volumes{'F'}=189.9;
	$volumes{'P'}=112.7;
	$volumes{'S'}=89;
	$volumes{'T'}=116.1;
	$volumes{'W'}=227.8;
	$volumes{'Y'}=193.6;
	$volumes{'V'}=140;

	#C-beta branching (number of non-hydrogen substituents attached to C-beta carbon)
	$cbetabranching{'A'}=1;
	$cbetabranching{'R'}=1;
	$cbetabranching{'D'}=1;
	$cbetabranching{'N'}=1;
	$cbetabranching{'C'}=1;
	$cbetabranching{'E'}=1;
	$cbetabranching{'Q'}=1;
	$cbetabranching{'G'}=1;
	$cbetabranching{'H'}=1;
	$cbetabranching{'I'}=2; #2xC
	$cbetabranching{'L'}=1;
	$cbetabranching{'K'}=1;
	$cbetabranching{'M'}=1;
	$cbetabranching{'F'}=1;
	$cbetabranching{'P'}=1;
	$cbetabranching{'S'}=1;
	$cbetabranching{'T'}=2; #1xC + 1xO
	$cbetabranching{'W'}=1;
	$cbetabranching{'Y'}=1;
	$cbetabranching{'V'}=2; #2xC

	#Kyte-Doolittle hydrophobicity
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

	$formalcharge{'D'}=-1;
	$formalcharge{'E'}=-1;
	$formalcharge{'H'}=1;
	$formalcharge{'K'}=1;
	$formalcharge{'R'}=1;
	$formalcharge{'I'}=0;
	$formalcharge{'V'}=0;
	$formalcharge{'L'}=0;
	$formalcharge{'F'}=0;
	$formalcharge{'C'}=0;
	$formalcharge{'M'}=0;
	$formalcharge{'A'}=0;
	$formalcharge{'G'}=0;
	$formalcharge{'T'}=0;
	$formalcharge{'S'}=0;
	$formalcharge{'W'}=0;
	$formalcharge{'Y'}=0;
	$formalcharge{'P'}=0;
	$formalcharge{'Q'}=0;
	$formalcharge{'N'}=0;
	}
