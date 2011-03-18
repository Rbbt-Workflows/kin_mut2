#tgi_bio..pm
#GENERAL PURPOUSE MODULE - PERMORFS BASIC BIOLOGICAL TASKS

#AUTHOR: TXEMA GONZALEZ IZARZUGAZA (2006)
#EMAIL: jmgonzalez@cnio.es
#STRUCTURAL BIOINFORMATICS UNIT - SPANISH NATIONAL CANCER RESEARCH CENTER (CNIO)

package tgi_bio;
require Exporter;
@ISA=qw(Exporter);
@EXPORT=qw(
	tgi_convert_aa_code_3to1
	tgi_convert_aa_code_1to3
	tgi_get_name_from_aa_code
	tgi_fasta2stringseq
	tgi_seq2array
	); 

use strict;

sub tgi_convert_aa_code_3to1{
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Converts a 3 letters aminoacid code to the corresponding 1 letter code.
	#INPUT: The 3 letters code
	#OUTPUT:The 1 letter code
	my $code3=lc($_[0]);
	my %codes;
	$codes{"ala"}="A";
	$codes{"arg"}="R";
	$codes{"asn"}= "N";
	$codes{"asp"}="D";
	$codes{"cys"}="C";
	$codes{"gln"}="Q";
	$codes{"glu"}="E";
	$codes{"his"}="H" ;
	$codes{"ile"}= "I";
	$codes{"gly"}= "G";
	$codes{"leu"}="L";
	$codes{"lys"}= "K";
	$codes{"met"}= "M";
	$codes{"phe"}= "F";
	$codes{"pro"}= "P";
	$codes{"ser"}= "S";
	$codes{"thr"}= "T";
	$codes{"trp"}= "W";
	$codes{"tyr"}="Y";
	$codes{"val"}= "V";
	return $codes{$code3};
	}


sub tgi_convert_aa_code_1to3{
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Converts a 1 letter aminoacid code to the corresponding 3 letters code.
	#INPUT: The 1 letter code
	#OUTPUT:The 3 letters code
	my $code1=uc($_[0]);
	my %codes;
	$codes{"A"}="ala";
	$codes{"R"}="arg";
	$codes{"N"}= "asn";
	$codes{"D"}="asp";
	$codes{"C"}="cys";
	$codes{"Q"}="gln";
	$codes{"E"}="glu";
	$codes{"H"}="his" ;
	$codes{"I"}= "ile";
	$codes{"G"}= "gly";
	$codes{"L"}="leu";
	$codes{"K"}= "lys";
	$codes{"M"}= "met";
	$codes{"F"}= "phe";
	$codes{"P"}= "pro";
	$codes{"S"}= "ser";
	$codes{"T"}= "thr";
	$codes{"W"}= "trp";
	$codes{"Y"}="tyr";
	$codes{"V"}= "val";
	return $codes{$code1};
	}
	
sub tgi_get_name_from_aa_code{
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Converts an aminoacid given as 1 or 3 letters code to its english name.
	#INPUT: The 1/3 letter code
	#OUTPUT:The english name
	my $code=uc($_[0]);
	my %codes;
	if ($code=~/^\s*\w\w\w\s*$/){
		$code=tgi_convert_aa_code_3to1($code);
		}
	$codes{"A"}="alanine";
	$codes{"R"}="arginine";
	$codes{"N"}= "asparagine";
	$codes{"D"}="aspartic acid";
	$codes{"C"}="cysteine";
	$codes{"Q"}="glutamine";
	$codes{"E"}="glutamic acid";
	$codes{"H"}="histidine" ;
	$codes{"I"}= "isoleucine";
	$codes{"G"}= "glycine";
	$codes{"L"}="leucine";
	$codes{"K"}= "lysine";
	$codes{"M"}= "metionine";
	$codes{"F"}= "phenylalanine";
	$codes{"P"}= "proline";
	$codes{"S"}= "serine";
	$codes{"T"}= "threonine";
	$codes{"W"}= "tryptophan";
	$codes{"Y"}="tyrosine";
	$codes{"V"}= "valine";
	return $codes{$code};
	}
	
sub tgi_fasta2stringseq {
	my $fasta=$_[0]; #Reference to the array with the fasta sequence
	my $seq;
	for (my $i=1;$i<scalar @$fasta;$i++) {
		$seq.=$$fasta[$i];
		}
	return $seq;
	}
	
sub tgi_seq2array {
	my $seq=$_[0]; #Sequence (protein or genomic)
	my @raw=split(/[\s]*/,$seq);
	my @seq;
	my $pos=0;
	foreach my $res (@raw) {
		$pos++; #So, first position will be 1.
		$seq[$pos]=$res;
		}
	return @seq;
	}
	
	
	
