#!/usr/bin/perl
###################################################################################################################
#########                                                                                                      ####
#########    DEVELOPED BY THIERES TAYRONI MARTINS DA SILVA AND LEONARDO VINICIUS DIAS DA SILVA                 ####
#########    contact thierestayroni@gmail.com leomarie@gmail.com                                               ####
#########                                                                                                      ####
###################################################################################################################
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;
use Bio::SeqIO::fastq;
use Cwd;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use POSIX qw(strftime);

my $fastq_inputs; # path to file contains all reads in fastq.gz
my $ref_seqs; # path to file contains the references sequences in fasta
my $len_cutoff; # min sequence length to be analysed
my $qua_cutoff; # min quality fastq sequence to be analysed
my $max_cpu; #numer of cpus to be used
my $translate; # translate reads to be analysed (yes or not)
my $codon_table; # only if translate is select
my $percent_quality_cutoff; # percentage of bases with value above the quality cutoff
my $sequence_read;
GetOptions (
  "in=s" => \$fastq_inputs, 
  "ref_seqs=s" => \$ref_seqs,
  "len_cutoff=i" => \$len_cutoff,
  "quality_cutoff=i" => \$qua_cutoff,
  "cpu=i"=> \$max_cpu,
  "codon_table=s" => \$codon_table,
  "translate" => \$translate,
  "quality_percent_cutoff=i" => \$percent_quality_cutoff,
  ""
) or die ("Error in command line arguments\n");

my $covid_path = dirname (rel2abs($0));
#print $covid_path;
##Path to binaries 
my $fastq_filter_path = "$covid_path/bin/fastq_quality_filter";
my $fastq2fasta_path = "$covid_path/bin/fastq_to_fasta";
my $mothur_path = "$covid_path/mothur/mothur";
my $fastq_join = "$covid_path/ea-utils/clipper/fastq-join";
my $diamond_path = "";
#

my $pm = Parallel::ForkManager->new($max_cpu);
my $date = strftime "%Y%m%d%H%M", localtime;
my $result_path = "$covid_path/results/$date";
mkdir "$result_path";

my @names;
my %file_name;
my $run_mode = "single";
opendir (DIR,$fastq_inputs);
my @files = readdir (DIR);
foreach my $file(@files){
  if ($file =~ m/_r2/i){
  $run_mode = "paired";
  }
}
foreach  my $file(@files){
  if (($file eq ".") || ($file eq "..")){
    next;
  }
  if ($file !~ /.fastq.gz$/){
    next;
  }
  if ($run_mode eq "paired"){
    my @splitname = split (/_/, $file);
    #print Dumper(@splitname);
    my $temp_name = "$splitname[0]_$splitname[1]_$splitname[2]";
    push (@names, $temp_name);
    mkdir "$result_path/$temp_name";
    $file_name{$file} = $temp_name;
    system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name\_$splitname[3].fastq";
  }
  else{
    my $temp_name = $file;
    $temp_name =~ s/.fastq.gz//g;
    $file_name{$file} = $temp_name;
    push (@names, $temp_name);
    mkdir "$result_path/$temp_name";
    system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name.fastq";
  }
}
$pm -> wait_all_children;

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

my @filter_names = uniq(@names);

if ($run_mode eq "paired"){
  JOIN:
  foreach my $name(@filter_names){
  my $pid = $pm -> start and next JOIN;
  #print $name."\n";
  system "$fastq_join $result_path/$name/$name\_R1.fastq $result_path/$name/$name\_R2.fastq -o $result_path/$name/$name\_ > /dev/null" ;
  my $seq_reverse = Bio::SeqIO->new ( -fortmat => 'fastq',
                                      -file => "$result_path/$name/$name\_un2");
			      print Dumper(\$seq_reverse);
  while (my $seq = $seq_reverse->next_seq()){
	  #print $seq."\n";
  }
  system "cat $result_path/$name/$name\_join $result_path/$name/$name\_un1 $result_path/$name/$name\_un2 > $result_path/$name/$name.fastq";
  $pm -> finish;
  }
}
$pm -> wait_all_children;
#READS:
#foreach my $file(@files){
#  my $pid = $pm -> start and next READS;
#  if (($file eq ".")|| ($file eq "..")){
#    $pm -> finish;
#    next;
#  }
#  my $path2file = "$fastq_inputs$file";
#  print "Descompactando arquivo $file\n";
#  my $fastq_name = $file;
#  $fastq_name =~ s/_R[1|2]_.+//;
#  system "gunzip -c $path2file > $result_path/$file_name{$file}/$fastq_name";
#  $pm -> finish;
#}
#$pm -> wait_all_children;

#print Dumper(@names);
QUALITY_FILTER:
foreach my $name(@names){
  my $pid = $pm -> start and next QUALITY_FILTER;
  #my $fastq_filter_path = "$covid_path/bin/fastq_quality_filter";
  my $seq_path = "$result_path/$name/$name.fastq";
  my $seq_out = "$result_path/$name/$name.fastq.filter";
  print "Filtrando arquivo $name pela qualidade do sequenciamento\n";
  system (" $fastq_filter_path -Q 33 -q $qua_cutoff -p $percent_quality_cutoff -i $seq_path -o $seq_out");
  $pm -> finish;
}
$pm -> wait_all_children;

FASTQ2FASTA:
foreach my $name(@names){
   my $pid = $pm -> start and next FASTQ2FASTA;
   my $seq_path = "$result_path/$name/$name.fastq.filter";
   my $seq_out = "$result_path/$name/$name.fastq.filter.fasta";
   print "Convertendo o arquivo $name de fastq para fasta\n";
   system ("$fastq2fasta_path -Q33 -i $seq_path -o $seq_out -n");
   $pm -> finish;
}
$pm -> wait_all_children;

UNIQUE:
foreach my $name(@names){
  my $pid = $pm -> start and next UNIQUE;
  my $filter_seq_path = "$result_path/$name/$name.fastq.filter.fasta";
  print "Sumarizando sequencias identicas do arquivo $name.\n";
  system ("$mothur_path \"#set.logfile(name=silent);unique.seqs(fasta=$filter_seq_path, format=count)\" ");
  $pm -> finish;
}
$pm -> wait_all_children;


print "Obrigado por escolher a PHT\n";
