#tgi_basic..pm
#GENERAL PURPOUSE MODULE

#AUTHOR: TXEMA GONZALEZ IZARZUGAZA (2006)
#EMAIL: jmgonzalez@cnio.es
#STRUCTURAL BIOINFORMATICS UNIT - SPANISH NATIONAL CANCER RESEARCH CENTER (CNIO)

package tgi_basic;
require Exporter;
@ISA=qw(Exporter);
@EXPORT=qw(
	tgi_open
	tgi_create_dir
	tgi_create_file
	tgi_check_input
	tgi_error
	tgi_timestamp
	tgi_log_base
	tgi_abs
	tgi_print_array
	tgi_print_hash
	tgi_get_filename
	tgi_get_path
	tgi_minimum
	tgi_minimum_couple
	tgi_maximum
	tgi_maximum_couple
	tgi_mean
	tgi_median
	tgi_clean_dir
	tgi_print_exec_header
	tgi_print_exec_footer
	tgi_random_tag
	tgi_delete_repeats
	tgi_delete_file
	tgi_clean_temporary_filesystem
	tgi_shuffle
	); 

sub tgi_shuffle {
	  my @a=\(@_);
	  my $n;
	  my $i=@_;
	  map {
	    $n = rand($i--);
	    (${$a[$n]}, $a[$n] = $a[$i])[0];
	  } @_;
	}

sub tgi_open {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Opens a file as read only.
	#INPUT: The file path
	#OUTPUT: An array with the file
	my $file=$_[0];
	open (TGI_OPEN,"<$file") or die("ERROR: Can't open: $file\n") ;
	chomp(my @file=<TGI_OPEN>);
	close(TGI_OPEN);
	return @file;
	}

sub tgi_create_dir{
	my $dir=$_[0];
	unless (-e $dir) {
		system("mkdir $dir");
		}
	return 0;
	}
	
sub tgi_create_file{
	my $file=$_[0];
	unless (-e $file) {
		system("touch $file");
		}
	return 0;
	}

sub tgi_timestamp {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Returns the date and hour, to make a timestamp of an analysis
	#INPUT: NONE
	#OUTPUT: the date as:  (Mon Nov 13 18:54:58 CET 2006)
	my $date=`date`;
	chomp($date);
	return $date;
	}

sub tgi_get_filename{
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Gets the name of a file (excluding its path)
	#INPUT: The file longname (path+name)
	#OUTPUT: The file name
	my $fullfilename=$_[0];
	my @tmp=split(/\//,$fullfilename);
	my $max=scalar @tmp-1;
	my $filename=$tmp[$max];
	return $filename;
	}

sub tgi_get_path{
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Gets the path of a file 
	#INPUT: The file longname (path+name)
	#OUTPUT: The file path
	my $filename=$_[0];
	if ($filename=~/([\w\/]*[\w\W]+\/)/){#DIR/FILE
		return $1;
		}
	my $file=tgi_get_filename($filename);
	$file=~s/\///gi;
	if ($filename eq "/$file" ){#ROOT
		return "/";
		}
	return "." #FILE IN THIS DIRECTORY, but not specified
	}

sub tgi_check_input {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Checks the input and displays an error message with the usage
	#INPUT: 1) The expected number of inputs
	#	2) The usage of the script
	#OUTPUT: None
	my $n=$_[0]; 
	my $usage=$_[1];
	if (scalar @ARGV<$n) {
		&tgi_error("Not enough parameters\nUSAGE: $usage\nForcing Exit\!\n\n");
		}
	return 1;
	}

sub tgi_error {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Displays an error message and finishes the script
	#INPUT: The Error Message.
	#OUTPUT: None
	my $msg=$_[0] ;
	die ("ERROR: $msg") ;
	return 1;
	}

sub tgi_print_array {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Prints an array, each value in a line
	#INPUT: Array to be printed.
	#OUTPUT: None
	my @array=@_;
	my $lin;
	foreach $lin (@array) {
		print "$lin\n";
		}
	return;
	}

sub tgi_print_hash {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Prints a hash, each key => value in a line
	#INPUT: Reference to the HASH to be printed.
	#OUTPUT: None
	my $hash=$_[0];
	my $lin;
	foreach $lin (keys %$hash) {
		print "$lin => $$hash{$lin}\n";
		}
	return;
	}
	
sub tgi_log_base {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Calculates the log of a number in a base different that 'e'
	#INPUT: The new base, and the value.
	#OUTPUT: The log (base $base) of $value
	my ($base,$value)=@_;
	return log($value)/log($base);
	}
	
sub tgi_abs {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Calculates the absolute value of a number
	#INPUT: The value.
	#OUTPUT: The absolute $value
	my $value=$_[0];
	if ($value<0){
		$value*=-1;
		}
	return $value;
	}

sub tgi_minimum {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Returns the minimum value, given a list
	#INPUT: A set of values
	#OUTPUT: The minimum of the values
	my @values=@_;
	my @sorted=sort {$a <=> $b} @values;
	return $sorted[0];
	}

sub tgi_minimum_couple {
	my $max=$_[0];
	my $min=$_[1];
	
	if ($min>$max) {
		$max=$_[1];
		$min=$_[0];
		}
	
	return $min;
	}

sub tgi_maximum {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Returns the maximum value, given a list
	#INPUT: A set of values
	#OUTPUT: The maximum of the values
	my @values=@_;
	my @sorted=reverse sort {$a <=> $b} @values;
	return $sorted[0];
	}

sub tgi_maximum_couple {
	my $max=$_[0];
	my $min=$_[1];
	
	if ($min>$max) {
		$max=$_[1];
		$min=$_[0];
		}
	return $max;
	}
	
sub tgi_mean {
	my @values=@_;
	my $mean;
	my $sum;
	my $n=scalar @values;
	foreach (@values) {
		$sum+=$_;
		}
	if ($n!=0) {
		$mean=$sum/$n;
		}
	else {
		$mean="INF";
		}
	return $mean;
	}	

sub tgi_median {
	my $median;
	my @values= sort {$a <=> $b} @_;
	
	my $n=scalar @values;
	my $r=$n%2;
	
	if ($r==0) {
		my $center=$n/2;
		$median=tgi_mean($values[$center],$values[$center-1]);
		}
		
	else {
		my $center=($n/2)-0.5;
		$median=$values[$center+0.5];
		}
	return $median;
	}

sub tgi_clean_dir {
	#AUTHOR: Txema Glez Izarzugaza
	#DESCRIPTION: Removes all files in a folder, one by one
	#INPUT: A path
	#OUTPUT: NONE
	my $path=$_[0];
	opendir (DIR,$path);
	my @dir=readdir(DIR);
	close (DIR);
	foreach my $file (@dir){
		if ($file!~/^\.+$/) {
			system("rm -rf $path/$file");
			}
		}
	return 0;
	}

sub tgi_print_exec_header{
	system("clear");
	my $pl=tgi_get_filename($0);
	print STDERR "#################################################################\n";
	print STDERR "##\n";
	print STDERR "##\t\t$pl\n##\n";
	print STDERR "##\tAuthor:\tJose M. Gonzalez-Izarzugaza Martinez\n";
	print STDERR "##\tEmail:\tjmgonzalez\@cnio.es\n";
	my $date=tgi_timestamp();
	print STDERR "##\tStart:\t$date\n##\n";
	print STDERR "#################################################################\n";
	print STDERR "\n\n\n";
	}
	
sub tgi_print_exec_footer{
	my $pl=tgi_get_filename($0);	
	my $date=tgi_timestamp();
	print STDERR "\n\n#\tExecution of $pl finished properly\n";
	print STDERR "#\tEnd:\t$date\n";
	}

sub tgi_random_tag {
	srand (time ^ $$);
	my $tag=int(rand(10000));
	return $tag;
	}

sub tgi_delete_repeats {
	my @with_repeats=@_;
	my @without_repeats;
	my %temporary;
	foreach my $element (@with_repeats) {
		$temporary{$element}=1;
		}
	my @without_repeats=sort keys (%temporary);
	return @without_repeats;
	}
	
sub tgi_delete_file {
	my $file=$_[0];
	if (-e $file) {
		system("rm -rf $file &> /tmp/null");
		}
	return 0;
	}
	
sub tgi_clean_temporary_filesystem {
	my @tmp_files=@_;
	foreach my $tmp (@tmp_files) {
		tgi_delete_file($tmp);
		}
	return 0;
	}