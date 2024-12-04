#!/usr/bin/perl

# 09/2023

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use Statistics::Basic qw(:all nofill);
use List::Util qw(sum);

#~ use Array::Utils qw(:all);
#~ use List::Compare;
#~ use List::MoreUtils qw(uniq);
#~ use Storable;
#~ use feature 'unicode_strings';
#~ use utf8;

#~ use Math::Round qw(:all);


my %opts = ();
getopts( 'i:l:s:c:o:h', \%opts ) or do_help();
if (defined ($opts{i})){
	#~ print "-i $opts{i}\n";
}
else{
	do_help();
	exit;
}
if (defined ($opts{l})){
	#~ print "-l $opts{l}\n";
}
else{
	do_help();
	exit;
}
if (defined ($opts{s})){
	#~ print "-s $opts{s}\n";
}
else{
	do_help();
	exit;
}
if (defined ($opts{c})){
	#~ print "-c $opts{c}\n";
}
else{
	do_help();
	exit;
}
if (defined ($opts{o})){
	#~ print "-l $opts{l}\n";
}
else{
	do_help();
	exit;
}
if ( $opts{h} ) {
	do_help();
}


sub do_help {
printf "$0";
	printf "
	-i	Input coverage and FPK count for each windows
	-l	Input calculated median FPK value
    -s	Step to consider around a target region to not have any bases overlap to the target region (non-overlapped)
    -c  Expected minimum FPK value
	-o	Output file with corrected coverage and FPK
	
	Option -h : show the help\n"
}

#~ perl apply_local_lambda.pl -i genome_win_count.tsv -l medianFPK -s step -o genome_win_count_lambda_corrected.tsv

my $input_file=$opts{i};
my $output_file=$opts{o};

## default for expected coverage (10) --> need to apply lambda to X windows around the target one until in theroy we reach the expected coverage
my $expected_coverage=$opts{c};

my $medFPK=$opts{l};
my $step=$opts{s};


# Calculate the theoretical number of regions you must consider to reach the expected coverage
my $lambda_tot=int(($expected_coverage/$medFPK) + 0.9999);
# Check if the result is even; if so, add 1 to center around the target region
$lambda_tot++ if $lambda_tot % 2 == 0;
# Calculate the lambda to apply to each side of the target region
my $lambda_half=int($lambda_tot/2);

## No header
my @header=("SAF","Chromosome","Start","End","Strand","Length","CovInput","CovChiP","FPKInput","FPKChiP");		## No header
my ($col_saf, $col_chr, $col_start, $col_end, $col_strand, $col_len, $col_covinput, $col_covchip, $col_fpkinput, $col_fpkchip);

$col_saf=0;
$col_chr=1;
$col_start=2;
$col_end=3;
$col_strand=4;
$col_len=5;
$col_covinput=6;
$col_covchip=7;
$col_fpkinput=8;
$col_fpkchip=9;

open( IN, $input_file ) or die("$! : Fail to open IN file $input_file line ".__LINE__."\n");

my %hinput_cov=();
my %hinput_fpk=();

my %hline=();

while ( defined( my $line = <IN> ) ) {
	chomp($line);
	my @tab=split("\t",$line);
	my $chr=$tab[$col_chr];
	my $start=$tab[$col_start];
	my $end=$tab[$col_end];
	my $covinput=$tab[$col_covinput];
	my $fpkinput=$tab[$col_fpkinput];
	
	$hinput_cov{$chr}{$start}=$covinput;
	$hinput_fpk{$chr}{$start}=$fpkinput;

	$hline{$chr}{$start}=$line;
}
close(IN);

my %hinput_cov_corrected=();
my %hinput_fpk_corrected=();

# my $lambda=($lambda_half*$step)+1;
my $lambda=($lambda_half*$step);

my $lambda_twice=($lambda*2)+1;

foreach my $chr (keys %hinput_cov){
	
	my $reg_count=keys %{$hinput_cov{$chr}};		##nb_region per chr (actual size, 1-based)
	if($reg_count>$lambda_twice){
		
		my @input_cov=();
		foreach my $start (sort {$a <=> $b} keys %{$hinput_cov{$chr}}){
			push(@input_cov,$hinput_cov{$chr}{$start});                        ## @input_cov is 0-based
		}
		my @input_fpk=();
		foreach my $start (sort {$a <=> $b} keys %{$hinput_fpk{$chr}}){
			push(@input_fpk,$hinput_fpk{$chr}{$start});                        ## @input_fpk is 0-based
		}
		
		my $i=1;
		my $end_of_chr=($reg_count-$lambda);
		my $end_of_chr_reg=($reg_count-$lambda_twice);
		
		foreach my $start (sort {$a <=> $b} keys %{$hinput_cov{$chr}}){
			if($i<=$lambda){							## lambda unbalanced toward 3' (start of chr)
				
				my $tot_cov=0;
				my $tot_fpk=0;
				for (my $j=0; $j<$lambda_twice; $j++){
					$tot_cov+=$input_cov[$j];
					$tot_fpk+=$input_fpk[$j];
				}
				my $avg_cov=($tot_cov/$lambda_twice);
				$hinput_cov_corrected{$chr}{$start}=$avg_cov;
				my $avg_fpk=($tot_fpk/$lambda_twice);
				$hinput_fpk_corrected{$chr}{$start}=$avg_fpk;
				
			}
			elsif($i>=$end_of_chr){					## lambda unbalanced toward 5' (end of chr)
				
				my $tot_cov=0;
				my $tot_fpk=0;
				for (my $j=$end_of_chr_reg; $j<=$#input_cov; $j++){
					$tot_cov+=$input_cov[$j];
					$tot_fpk+=$input_fpk[$j];
				}
				my $avg_cov=($tot_cov/$lambda_twice);
				$hinput_cov_corrected{$chr}{$start}=$avg_cov;
				my $avg_fpk=($tot_fpk/$lambda_twice);
				$hinput_fpk_corrected{$chr}{$start}=$avg_fpk;
				
			}
			else{
				my $tot_cov=0;
				my $tot_fpk=0;									## lambda balanced
				
				#	$lambda (keys %hinput_cov) is 1-based and @hinput_cov ($j) is 0-based 
				my $start_reg=($i-$lambda)-1;
				my $end_reg=($i+$lambda)-1;

				for (my $j=$start_reg; $j<=$end_reg; $j++){
					$tot_cov+=$input_cov[$j];
					$tot_fpk+=$input_fpk[$j];
				}
				my $avg_cov=($tot_cov/$lambda_twice);
				$hinput_cov_corrected{$chr}{$start}=$avg_cov;
				my $avg_fpk=($tot_fpk/$lambda_twice);
				$hinput_fpk_corrected{$chr}{$start}=$avg_fpk;
			}
			$i++;
		}
	}
	else{											##Average for all if less than twice the lambda
		my @input_cov=values %{$hinput_cov{$chr}};
		my $avg_cov=(sum(@input_cov)/($#input_cov+1));
		
		my @input_fpk=values %{$hinput_fpk{$chr}};
		my $avg_fpk=(sum(@input_fpk)/($#input_fpk+1));
		
		foreach my $start (keys %{$hinput_cov{$chr}}){
			$hinput_cov_corrected{$chr}{$start}=$avg_cov;
			$hinput_fpk_corrected{$chr}{$start}=$avg_fpk;
		}
	}
}


open( OUT, ">".$output_file ) or die("$! : Fail to open OUT file $output_file line ".__LINE__."\n");
print OUT join("\t",@header)."\tCovInputCorrected\tFPKInputCorrected\n";
foreach my $chr (sort {$a cmp $b} keys %hinput_cov){
	foreach my $start (sort {$a <=> $b} keys %{$hinput_cov{$chr}}){
		print OUT $hline{$chr}{$start}."\t".$hinput_cov_corrected{$chr}{$start}."\t".$hinput_fpk_corrected{$chr}{$start}."\n";
	}
}
close(OUT);