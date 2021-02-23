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
use Bio::SeqIO;
use Bio::Seq::Quality;
use Cwd;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use POSIX qw(strftime);

my $fastq_inputs; # path to directory contains all files in fastq.gz
my $primer_input; #path to file contains all primers
my $ref_seqs; # path to file contains the references sequences in fasta
my $len_cutoff; # min sequence length to be analysed
my $qua_cutoff; # min quality fastq sequence to be analysed
my $weigth_cutoff; #min number of identical sequence to be analysed
my $max_cpu; #numer of cpus to be used
my $translate; # translate reads to be analysed (yes or not)
my $codon_table; # only if translate is select
my $percent_quality_cutoff; # percentage of bases with value above the quality cutoff
my $sequence_read;
my $depara; #Illumina Id name to real name;
GetOptions (
  "in=s" => \$fastq_inputs, 
  "primers=s" => \$primer_input,
  "ref_seqs=s" => \$ref_seqs,
  "len_cutoff=i" => \$len_cutoff,
  "quality_cutoff=i" => \$qua_cutoff,
  "cpu=i"=> \$max_cpu,
  "codon_table=s" => \$codon_table,
  "translate" => \$translate,
  "quality_percent_cutoff=i" => \$percent_quality_cutoff,
  "weigth=i" => \$weigth_cutoff,
  "depara=s" => \$depara,
) or die ("Error in command line arguments\n");

my $covid_path = dirname (rel2abs($0));
#print $covid_path;
##Path to binaries 
my $fastq_filter_path = "$covid_path/bin/fastq_quality_filter";
my $fastq2fasta_path = "$covid_path/bin/fastq_to_fasta";
my $mothur_path = "$covid_path/mothur/mothur";
my $fastq_join = "$covid_path/ea-utils/clipper/fastq-join";
my $diamond_path = "$covid_path/diamond/bin/diamond";
my $blast_makedb_path = "$covid_path/ncbi-blast-2.11.0+/bin/makeblastdb";
my $blastn_path = "$covid_path/ncbi-blast-2.11.0+/bin/blastn";
my $blast_formatter = "$covid_path/ncbi-blast-2.11.0+/bin/blast_formatter";
#

my $pm = Parallel::ForkManager->new($max_cpu);
my $date = strftime "%Y%m%d%H%M", localtime;
my $result_path = "$covid_path/results/$date";
mkdir "$result_path";
mkdir "$result_path/final_resuts";

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
    print "Extraindo arquivo $temp_name\n";
    system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name\_$splitname[3].fastq";

  }
  else{
    my $temp_name = $file;
    $temp_name =~ s/\_R1.*\.fastq\.gz//i;
    $file_name{$file} = $temp_name;
    push (@names, $temp_name);
    mkdir "$result_path/$temp_name";
    print "Extraindo arquivo $temp_name\n";
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
  $pm -> finish;
  }
  $pm -> wait_all_children;
  
  foreach my $name(@filter_names){
    open (IN, "$result_path/$name/$name\_un2");
    open (OUT, ">>/$result_path/$name/$name\_un2rev");
    my $i = 1;
    while (my $line = <IN>) {
      chomp $line;
      if ($i == 1){
	my @fields = split (/\s/, $line);
	#if ($line !~ m/^@/);
        print OUT "$fields[0]:2 $fields[1]\n";
        $i++;
      }
      elsif ($i == 2){
        print OUT scalar reverse ($line);
        print OUT"\n";
        $i++;
      }
      elsif ($i == 3){
        print OUT "$line \n";
        $i++;
      }
      elsif ($i == 4){
        print OUT scalar reverse ($line);
        print OUT "\n";
      
        $i = 1;
      }
    }
  close (IN);
  close (OUT);
  system "cat $result_path/$name/$name\_join $result_path/$name/$name\_un1 $result_path/$name/$name\_un2rev > $result_path/$name/$name.fastq";
  } 
}


QUALITY_FILTER:
foreach my $name(@filter_names){
  my $pid = $pm -> start and next QUALITY_FILTER;
  #my $fastq_filter_path = "$covid_path/bin/fastq_quality_filter";
  my $seq_path = "$result_path/$name/$name.fastq";
  if (!-s $seq_path){
    $pm->finish;
    next;
  }
  my $seq_out = "$result_path/$name/$name.qfilter.fastq";
  print "Filtrando arquivo $name pela qualidade do sequenciamento\n";
  system (" $fastq_filter_path -Q 33 -q $qua_cutoff -p $percent_quality_cutoff -i $seq_path -o $seq_out");
  $pm -> finish;
}
$pm -> wait_all_children;

FASTQ2FASTA:
foreach my $name(@filter_names){
   my $pid = $pm -> start and next FASTQ2FASTA;
   my $seq_path = "$result_path/$name/$name.qfilter.fastq";
   if (!-s $seq_path){
    $pm->finish;
    next;
   }
   my $seq_out = "$result_path/$name/$name.qfilter.fasta";
   print "Convertendo o arquivo $name de fastq para fasta\n";
   system ("$fastq2fasta_path -Q33 -i $seq_path -o $seq_out -n");
   $pm -> finish;
}
$pm -> wait_all_children;


my %sequencew;
my %primerfr;
my @primers;
open (PRIMERS, "$primer_input") or die ("Não consegui abrir o arquivo de primer $primer_input\n");
while (<PRIMERS>){
  my $line = $_;
  chomp $line;
  
  my @fields = split (/\s/, $line);
  if ($fields[0] eq "R"){
    $fields[1] = reverse scalar $fields[1];
  }
  push (@primers, $fields[1]);
  $primerfr{$fields[1]} = $fields[0];
}
close (PRIMERS);
$/=">";

LENGTH:
foreach my $name(@filter_names){
  my $pid = $pm -> start and next LENGTH;
  if (!-s "$result_path/$name/$name.qfilter.fasta"){
    $pm-> finish;
    next;
  }
  open (IN, "$result_path/$name/$name.qfilter.fasta");
  open (OUT, ">>$result_path/$name/$name.qfilter.pfilter.fasta");
  while (<IN>){
    my $line= $_;
    
    my @fasta_line = split (/\n/, $line);
    my $id = $fasta_line[0];
    my $seq = $fasta_line[1];
    if (!$seq){
      next;
    }
    
    my $seq_length = length($seq);
    if ($seq_length < $len_cutoff){
      next;
    }
    
    foreach my $primer(@primers){
      if ($seq =~ /$primer/){
        if ($primerfr{$primer} eq "F"){
          $seq =~ s/^.*$primer//;
        }
	      else{
	        $seq =~ s/$primer.*$//;
	      }
      }
    }
    print OUT ">$id\n$seq\n";
  }
  close (IN);
  close (OUT);
  $pm->finish;
}
$pm -> wait_all_children;

$/="\n";

UNIQUE:
foreach my $name(@filter_names){
  my $pid = $pm -> start and next UNIQUE;
  my $filter_seq_path = "$result_path/$name/$name.qfilter.pfilter.fasta";
  if (!-s $filter_seq_path){
    $pm->finish;
    next;
  }
  print "Sumarizando sequencias identicas do arquivo $name\n";
  system ("$mothur_path \"#set.logfile(name=silent);unique.seqs(fasta=$filter_seq_path, format=count)\" > /dev/null");
  $pm -> finish;
}
$pm -> wait_all_children;

foreach my $name(@filter_names){
  my $file_count_path = "$result_path/$name/$name.qfilter.pfilter.count_table";
  if (!-s $file_count_path){
    next;
  }
  open (COUNT, "$file_count_path");
  while (<COUNT>){
    my $line = $_;
    chomp $line;
    my @fields = split (/\s/, $line);
    $sequencew{$name}{$fields[0]} = $fields[1];
  }
  close (COUNT);
}
#print Dumper(\%sequencew);
#print "\n$weigth_cutoff";
$/=">";

WEIGTH:
foreach my $name(@filter_names){
  my $pid = $pm->start and next WEIGTH;
  my $fasta_file = "$result_path/$name/$name.qfilter.pfilter.unique.fasta";
  if (!-s $fasta_file){
    $pm->finish;
    next;
  }
  my $fasta_out = "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta";

  open (IN, $fasta_file) or die;
  open (OUT, ">>$fasta_out");
  while (<IN>){
    my $line = $_;
    my @field_line = split (/\n/, $line);
    my $id = $field_line[0];
    my $seq = $field_line[1];
    my @id_fields = split (/\s/, $id);
    #print "o ID é $id\n";
    #print ("$sequencew{$name}{$id} heheu\n");
    #print $id."\t lalala\n";
    #print $seq."\t euheuhe\n";
    if (!$seq){
      next;
    }
    if ($sequencew{$name}{$id_fields[0]} < $weigth_cutoff){
      next;
    }
    print OUT (">$id_fields[0]\t$sequencew{$name}{$id_fields[0]}\n$seq\n");
  }
  close (IN);
  close (OUT);
  $pm -> finish;
}
$pm-> wait_all_children;
$/="\n";
BLAST:
foreach my $name(@filter_names){
  my $pid = $pm -> start and next BLAST;
  my $fasta = "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta";
  my $database = "$result_path/$name/$name.db.ndb";
  if (!-s $fasta){
    $pm->finish;
    next;
  }
  print "Rodando blast pro arquivo $name\n";
  system "$blast_makedb_path -in $fasta -dbtype nucl -title $name.db -out $database > /dev/null";
  system "$blastn_path -task blastn -query $ref_seqs -evalue 1 -db $database -outfmt 11 -out $result_path/$name/$name\_blast_out";
  system "$blast_formatter -archive $result_path/$name/$name\_blast_out -outfmt \"6 qseqid sseqid pident qcovs evalue\" -out $result_path/$name/$name\_blast.tsv";
  system "$blast_formatter -archive $result_path/$name/$name\_blast_out -outfmt 2 -out $result_path/$name/$name\_blast_alg";

  $pm -> finish;
}
$pm-> wait_all_children;

my @ids_ref;
my $in = Bio::SeqIO->new (-file =>"$ref_seqs", -format => "Fasta");
while (my $seq_ref = $in->next_seq() ){
  my $id = $seq_ref->id();
  push (@ids_ref, $id);
}
my $header = "Indivíduo";
foreach my $id_ref(@ids_ref){
  $header = "$header,\"$id_ref\nResultado/Identidade/Cobertura/Peso\"";
}
open (DEPARA, $depara) or die;
my %realnames;
while (<DEPARA>){
  chomp $_;
  my @fields = split (/,/, $_);
  $realnames{$fields[0]} = $fields[1];
}
close (DEPARA);
open (OUTCSV, ">>$result_path/final_resuts/resultado_final.csv") or die ("Não consegui crar  o arquivo resultado_final.csv\n");
print OUTCSV "$header\n";
my %results_ref;
my %status;
foreach my $name(@filter_names){
  print OUTCSV "$realnames{$name}";
  my $idvalue;
  my $idcount =0;
  my $line2print;
  if (-e "$result_path/$name/$name\_blast.tsv"){
    foreach my $id(@ids_ref){
      open (BLASTTSV, "$result_path/$name/$name\_blast.tsv"); #or die ("Não consegui abri o arquivo $result_path/$name/$name\_blast.tsv");
      $status{$name}{$id} = 0; # 0 = negativo;
      foreach my $line(<BLASTTSV>){
        if ($line =~ /^$id/){
          my @fields = split (/\t/, $line);
          if ($fields[2] > 90){
            $status{$name}{$id} = 1; #1 = positivo;
            if ($idcount<1){
              $idvalue .= "+/$fields[2]\%/$fields[3]\%/$sequencew{$name}{$fields[1]}";
              $idcount++;
            }
            else{
              $idvalue.= "\n+/$fields[2]\%/$fields[3]\%/$sequencew{$name}{$fields[1]}";
            }
          }
        }else{next}
      
      }
      if ($status{$name}{$id} == 0){
        $idvalue = "-";
      }
      if ($idvalue){
          $line2print.= ",\"$idvalue\"";
        }
      $idcount=0;
      $idvalue="";
      close(BLASTTSV);
    }
  } else{
    $line2print.=",indefinido";
  }
  $line2print.= "\n";
  print OUTCSV "$line2print";
  $line2print="";
}

close (OUTCSV);

ALG:
foreach my $name(@filter_names){
  my $pid = $pm-> start and next ALG;
  my $blast_alg = "$result_path/$name/$name\_blast_alg";
  my %sequencescout;
  if (!-s $blast_alg){
    print "LALALA\n";
    $pm -> finish;
    next;
  }
  open (IN, $blast_alg) or die("Não consegui abrir o arquivo $blast_alg\n");
  open (FASTA, "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta") or die;
  open (OUT, ">>$result_path/final_resuts/$name.alg");
  my $i = 0;
  while (<FASTA>){
    chomp $_;
    print "$_\n";
      $sequencescout{$i} = $_;
  }
  print Dumper ("\%sequencescout\n");
  while (<IN>){
    chomp $_;
  }
}
print "Obrigado por escolher a PHT\n";