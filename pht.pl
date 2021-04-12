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
my $percent_quality_cutoff; # percentage of bases with value above the quality cutoff
my $depara; #Illumina Id name to real name;
my $ident_cutoff; # alignment identity percentage;
my $cov; # alignment coverage percentage;
GetOptions (
  "in=s" => \$fastq_inputs, 
  "primers=s" => \$primer_input,
  "ref_seqs=s" => \$ref_seqs,
  "len_cutoff=i" => \$len_cutoff,
  "quality_cutoff=i" => \$qua_cutoff,
  "cpu=i"=> \$max_cpu,
  "quality_percent_cutoff=i" => \$percent_quality_cutoff,
  "weigth=i" => \$weigth_cutoff,
  "depara=s" => \$depara,
  "ident_cutoff=i" => \$ident_cutoff,
  "cov=i"=> \$cov,
) or die ("Error in command line arguments\n");

my $covid_path = dirname (rel2abs($0));
#print $covid_path;
##Path to binaries 
my $fastq_filter_path = "$covid_path/dependencias/bin/fastq_quality_filter";
my $fastq2fasta_path = "$covid_path/dependencias/bin/fastq_to_fasta";
my $mothur_path = "$covid_path/dependencias/mothur/mothur";
my $fastq_join = "$covid_path/dependencias/ea-utils/clipper/fastq-join";
#my $diamond_path = "$covid_path/diamond/bin/diamond";
my $blast_makedb_path = "$covid_path/dependencias/ncbi-blast-2.11.0+/bin/makeblastdb";
my $blastn_path = "$covid_path/dependencias/ncbi-blast-2.11.0+/bin/blastn";
my $blast_formatter = "$covid_path/dependencias/ncbi-blast-2.11.0+/bin/blast_formatter";
#

my $pm = Parallel::ForkManager->new($max_cpu);
my $date = strftime "%Y%m%d%H%M", localtime;
print "O id da execução é $date\n";
my $result_path = "$covid_path/resultados/$date";
mkdir "$result_path";

mkdir "$result_path/resultados_finais";


my @names;
my %file_name;
my $run_mode = "single";
my $flag = 0;
opendir (DIR,$fastq_inputs);
my @files = readdir (DIR);
foreach my $file(@files){
  if ($file =~ m/_r2/i){
  $run_mode = "paired";
  }
  if ($file !~ m/_r\d/){
    $flag = 1;
  } 
}
foreach  my $file(@files){
  if (($file eq ".") || ($file eq "..")){
    next;
  }
  if ($file !~ /.fastq.gz$/){
    next;
  }
  if ($flag == 0){
    if ($run_mode eq "paired"){
      my @splitname = split (/_/, $file);
      #print Dumper(@splitname);
      my $i=0;
      my $r;
      foreach my $split(@splitname){
        if ($split =~ /r\d/i){
          $r = $i;
        }
        $i++;
      }
      my @tempslipt = split (/_R/, $file);
      my $temp_name = $tempslipt[0];
      push (@names, $temp_name);
      mkdir "$result_path/$temp_name";
      $file_name{$file} = $temp_name;
      print "Extraindo arquivo $file\n";
      system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name\_$splitname[$r].fastq";

    }
    else{
      my $temp_name = $file;
      $temp_name =~ s/\_R1.*\.fastq\.gz//i;
      $file_name{$file} = $temp_name;
      push (@names, $temp_name);
      mkdir "$result_path/$temp_name";
      print "Extraindo arquivo $file\n";
      system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name.fastq";
    }
  }
  if ($flag == 1){
    my $temp_name = $file;
    $temp_name =~ s/\.fastq\.gz//;
    $file_name{$file} = $temp_name;
    push (@names, $temp_name);
    mkdir "$result_path\/$temp_name";
    print "Extraindo arquivo $file\n";
    system "gunzip -c $fastq_inputs$file > $result_path/$temp_name/$temp_name.fastq";
  }
}
$pm -> wait_all_children;
close (DIR);

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
opendir (PRIMERS_DIR, "$primer_input") or die ("Não consegui abrir o arquivo de primer $primer_input\n");
my @primers_files = readdir(PRIMERS_DIR);
foreach  my $p_file(@primers_files){
  #print "$p_file\n";
  if (($p_file eq ".") || ($p_file eq "..")){
    next;
  }
  if ($p_file =~ /.csv/){
    open (PRIMERS, "$primer_input/$p_file") or die ("Não cnosegui abir o arquivo de primer $primer_input/$p_file\n");
    while (<PRIMERS>){
      my $line = $_;
      chomp $line;
  
      my @fields = split (',', $line);
      if ($fields[0] eq "R"){
        $fields[1] = reverse scalar $fields[1];
      }
      push (@primers, $fields[1]);
      $primerfr{$fields[1]} = $fields[0];
    }
  }
  close (PRIMERS);
}
close (PRIMERS_DIR);
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

#get fasta references file 
my $ref_path;
opendir (REFS, $ref_seqs) or die; 
my @ref_files = readdir (REFS);
foreach my $ref_file(@ref_files){
  if (($ref_file eq ".") || ($ref_file eq "..")){
    next;
  }
  if (($ref_file =~ /.fasta/) || ($ref_file =~ /.fa/)){
  $ref_path = "$ref_seqs/$ref_file";
  }
}


BLAST:
foreach my $name(@filter_names){
  my $pid = $pm -> start and next BLAST;
  my $fasta = "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta";
  my $database = "$result_path/$name/$name.db.ndb";
  if (!-s $fasta){
    $pm->finish;
    next;
  }
  print "Executando blast em $name\n";
  system "$blast_makedb_path -in $fasta -dbtype nucl -title $name.db -out $database > /dev/null";
  system "$blastn_path -task blastn -query $ref_path -evalue 1 -db $database -outfmt 11 -out $result_path/$name/$name\_blast_out";
  system "$blast_formatter -archive $result_path/$name/$name\_blast_out -outfmt \"6 qseqid sseqid pident qcovs evalue\" -out $result_path/$name/$name\_blast.tsv";
  system "$blast_formatter -archive $result_path/$name/$name\_blast_out -outfmt 2 -out $result_path/$name/$name\_blast_alg";

  $pm -> finish;
}
$pm-> wait_all_children;

my @ids_ref;
my %querys;
my $q = "Query_";
my $count_query = 1;
my $in = Bio::SeqIO->new (-file =>"$ref_path", -format => "Fasta");
while (my $seq_ref = $in->next_seq() ){
  my $id = $seq_ref->id();
  push (@ids_ref, $id);
}
foreach my$ref(@ids_ref){
  $querys{"$q$count_query"} = $ref;
  $count_query++;
}
my $header = "Indivíduo";
foreach my $id_ref(@ids_ref){
  $header = "$header,\"$id_ref\nResultado/Identidade/Cobertura/Peso/SeqID\"";
}

opendir (DEPARA_DIR, $depara) or die;

my %realnames;
my @d_files = readdir(DEPARA_DIR);
foreach my $d_file(@d_files){
  if (($d_file eq ".") || ($d_file eq "..")){
    next;
  }
  if ($d_file =~ /.csv/){
    system ("dos2unix -q $depara$d_file");
    open (DEPARA, "$depara/$d_file") or die ("Não cnosegui abir o arquivo de primer $depara/$d_file\n");
  }

  while (<DEPARA>){
    chomp $_;
    my @fields = split (',', $_);
    $realnames{$fields[0]} = $fields[1];
  }
  
}
close (DEPARA);
close (DEPARA_DIR);
open (OUTCSV, ">>$result_path/resultados_finais/resultado_final.csv") or die ("Não consegui crar  o arquivo resultado_final.csv\n");
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
          if (($fields[2] >= $ident_cutoff)&&($fields[3] >= $cov)){
            $status{$name}{$id} = 1; #1 = positivo;
            if ($idcount<1){
              $idvalue .= "+ / $fields[2]\% / $fields[3]\% / $sequencew{$name}{$fields[1]} / $fields[1]";
              $idcount++;
            }
            else{
              $idvalue.= "\n+ / $fields[2]\% / $fields[3]\% / $sequencew{$name}{$fields[1]} / $fields[1]";
            }
          }
        
          if (($fields[2] < $ident_cutoff)||($fields[3] < $cov)){
            if ($idcount<1){
              $idvalue .= "- / $fields[2]\% / $fields[3]\% / $sequencew{$name}{$fields[1]} / $fields[1]";
              $idcount++;
            }
            else{
              #print "Estou aqui\n";
              $idvalue.= "\n- / $fields[2]\% / $fields[3]\% / $sequencew{$name}{$fields[1]} / $fields[1]";
            }
          }
        }
      }
      #if ($status{$name}{$id} == 0){
      #  $idvalue = "-";
      #}
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

system "unix2dos $result_path/resultados_finais/resultado_final.csv";

my %sequencescout;
foreach my $name(@filter_names){
  if (!-s "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta"){
    next;
  }
  open (FASTA, "$result_path/$name/$name.qfilter.pfilter.unique.final.fasta") or die;
  my $i = 0;
  while (<FASTA>){
    my $line = $_;
    if ($line =~ /^>/){
      #print $line."\n";
      $sequencescout{$name}{$i} = $line;
      $i++;
    }
  }
}


ALG:
foreach my $name(@filter_names){
  my $pid = $pm-> start and next ALG;
  my $blast_alg = "$result_path/$name/$name\_blast_alg";
  if (!-s $blast_alg){
    #print "LALALA\n";
    $pm -> finish;
    next;
  }
  open (IN, $blast_alg) or die("Não consegui abrir o arquivo $blast_alg\n");
  open (OUT, ">>$result_path/resultados_finais/$name.txt");
  print OUT "Inicio\tSequência\tFim\tID\n";
  my $i =0;
  while (<IN>){
    chomp $_;
    if ($_ =~ /^Query\_/ ){
      if ($i>0){
        print OUT "\n";
      }
      my @fields = split (/\s/,$_);
      my $id = shift @fields;
      my $id_convert = %querys{$id};
      push (@fields, $id_convert);
      my $ref = join (' ', @fields);
      $ref =~ s/^\s*//;
      print OUT "$ref\n";
      $i++;
    }
    if ($_ =~ /^\d/ ){
      my @fields = split (/\s/,$_);
      my $id = $sequencescout{$name}{$fields[0]};
      chomp $id;
      $id =~ s/^>//;
      $id =~ s/\t\d*$//;
      shift @fields;
      push (@fields, $id);
      my $seq = join (' ', @fields);
      $seq =~ s/^\s*//;
      print OUT "$seq\n";
    }

    #if ($_ =~ /^Query\_/ ){
     # if ($i>0){
     #   print OUT "\n";
     # }
     # $_ =~ s/Query/ref/;
     # my @fields = split (/\s/,$_);
     # my $info = shift @fields;
     # push (@fields, $info);
     # my $ref;
     # foreach my $field(@fields){
     #   if ($field){
     #     $ref.="$field ";
     #   }
     # }
     # $ref =~ s/\s/\t/g;
     # print OUT "$ref\n";
     # $i++;
    #}
    #if ($_ =~ /^\d/){
    #  my @fields = split (/\s/, $_);
    #  #print Dumper (\@fields);
    #  my $id = $sequencescout{$name}{$fields[0]};
    #  chomp $id;
    #  $id =~ s/^>//;
    #  $id =~ s/\t\d*$//;
    #  $fields[0] = $id;
    #  my $info = shift @fields;
    #  push (@fields, $info);
    # # print Dumper (\@fields);
    #  my $seq;
    #  foreach my $field(@fields){
    #    if ($field){
    #      $seq.="$field ";
    #    }
    #  }
    #  $seq =~ s/\s/\t/g;
    #  print OUT "$seq\n";
    #} 
  }
  $pm->finish;
}
$pm -> wait_all_children;

print "Obrigado por escolher a PHT\n";
