#!/usr/bin/perl -w
# encoding: utf-8

use strict;
use warnings;

# simuCSVtoFASTA.pl
# A simple perl script to change populations experted in csv files exported from simuPOP, to
# a set of sequences in FASTA format.
# Based on script by Mark McMullan and modified by Ben Ward.

# This script assumes populations will be diploid.

# Takes simuPOP output, sequences in CSV format and converts them into fasta.
# The command used to get sequences within simupop is:
# sim.utils.Exporter(format='csv', output='outputfile.csv')

# Command Line ARGS:
# 0). Filename.csv, 
# 1). Name of sequences, 
# 2). Number of sequences sampled.


# We get the file to use from the command line argument, close if can't open it.
print STDERR "PARSING $ARGV[0]\n";
unless (open(CSV, $ARGV[0])){
  print "\nCan't open the .csv file\nClosing script\n";
  exit;
}
# Modify the output file name to a fasta name.
my $filename = $ARGV[0];
my $i;
for ($i = 0; $i < 4; $i ++){
  chop $filename;
}
my $output = "$filename.fas";
my $one = 'PloidOne';
my $two = 'PloidTwo';
my $sep1file = "$filename$one.fas";
my $sep2file = "$filename$two.fas";
# Open the output file.
open(OUTFILE, ">$output");
print OUTFILE "";
close(OUTFILE);
# Get the names of the sequences from the command line.
my $arrow = '>';
my $sequenceLabels = $ARGV[1];
print "The name to give to the sequences is: $sequenceLabels\n";
my $undscr ='_';
my $seq = "$sequenceLabels$undscr";
my $Individual = '_Ind_';
my $Number = 0;
my $Population = 'P';
my $a = 'a'; # This script assumes diploid populations.
my $b = 'b'; # This script assumes diploid populations.
my $counter = -1; # Clicks pop number ($num) up by 1
my $num = 1;
my $line = <CSV>; # Dispose of header line
# Get the number of sequences sampled from the command line.
my $nummax = $ARGV[2];
print "The number of sequences to sample is: $nummax\n";
#Load subsequent lines into array after removing whitespace.
foreach $line (<CSV>){
  $line =~ s/\s//g;
  $Number = $Number+1;
  $counter = $counter+1;
  if ($counter > $nummax){
    $counter = 1;
  }
  if ($counter == $nummax){
    $num = $num+1;
  }
  my @fasta12 = split(',', $line);

  # Get rid of the sex and effect fields in the sequence.
  shift @fasta12;
  shift @fasta12;

  # Split fasta12 into fasta1 (odd array entries) and fasta2 (even array entries).
  my $nuc;				
  my @fasta1;
  my @fasta2;
  my $length = scalar(@fasta12);
  print $length;
  # Divy out the base pairs / alleles into fasta1 or fasta2 i.e. choromosome set 1 or set 2.
  for ($i = 0; $i<$length; $i++){	# Check that I don't need to add 1 to $length
    my $oddevn = $i % 2;		# Remainder 0=even (or zero) remainder 1=odd
    $nuc = $fasta12[$i];
    if ($oddevn == 0){
      push(@fasta1, $nuc);
    }
    else{
      push(@fasta2, $nuc);
    }
  }
	
  # Now convert the arrays to strings.
  my $fasta1str = join('',@fasta1);
  my $fasta2str = join('',@fasta2);
  $fasta1str =~ s/0/A/g;
  $fasta1str =~ s/1/T/g;
  $fasta1str =~ s/2/C/g;
  $fasta1str =~ s/3/G/g;
  $fasta2str =~ s/0/A/g;
  $fasta2str =~ s/1/T/g;
  $fasta2str =~ s/2/C/g;
  $fasta2str =~ s/3/G/g;
  my $fastaName1 = "$arrow$seq$Population$num$Individual$Number$a";
  my $fastaName2 = "$arrow$seq$Population$num$Individual$Number$b";
  open(OUTFILE, ">>$output");
  print OUTFILE "$fastaName1\n$fasta1str\n$fastaName2\n$fasta2str\n";
  close(OUTFILE);
  open(PLOID1, ">>$sep1file");
  print PLOID1 "$fastaName1\n$fasta1str\n";
  close(PLOID1);
  open(PLOID2, ">>$sep2file");
  print PLOID2 "$fastaName2\n$fasta2str\n";
}


















