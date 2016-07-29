#!/usr/bin/perl -w
#
# 

# v3 2016.5.31 switch the case control position in list.txt
# Last modified by: Lance Chang (2016.5.23)
# Last modified by: Lance Chang (2016.1.7)
#
# Copyright belong to Chen Yuan Liu 
# Author: Chen-Yuan Liu
# Date: 10-28, 2015
# The Core for this pipeline is RSEM
#  
# Function:
#
#

$main = new RNASEQCore;
$main-> parse_command_line;
$main-> CalculateExpression;
$main-> GenerateDataMatrix;
$main-> ebseq;
# $main-> DEGexp;


package RNASEQCore;
use strict;
use warnings;
use Getopt::Long;
use threads;
# use vars qw($app %parameters @ctrldup @casedup @controls @cases $out_CE $out_DM $out_DEGseq $out_ebseq $ngvec $reference);

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  $self->init;
  return $self;
}
sub init {
  my $self = shift;
  our $ini;
  our $lst;
  our $parameters;
  
  GetOptions(
    'lst=s'    => \$lst,
    'ini=s'      => \$ini,
  );
  
  $parameters->{list}   = $lst;
  $parameters->{ini}    = $ini;   
  # print scalar(@ARGV)."\n";
  
  unless(-e $ini){    print usage();    exit(1);  } 

  our $output = "output";
  our $out_CE = "$output/CalculateExpression";
  our $out_DM = "$output/GenerateDataMatrix";  
  our $out_ebseq = "$output/ebseq";
  # our $out_DEGseq  = "$output/DEGseq";
 
  unless (-e $out_CE) {system ("mkdir -p $out_CE");}
  unless (-e $out_DM) {system ("mkdir -p $out_DM");}
  unless (-e $out_ebseq) {system ("mkdir -p $out_ebseq");}
  # unless (-e $out_DEGseq) {system ("mkdir -p $out_DEGseq");}  
}


sub parse_command_line {
  my $self = shift;
  our $parameters;
  our @sets;
  our $samples; # reference of %samples
  
  my $ln;


  my $lst=$parameters->{list};
  my $ini=$parameters->{ini};  

  open (IN, "<$ini") or die "Can't read $ini file";  
  while (<IN>) {
    $ln = $_;
    chomp($ln);
	
    if ($ln=~/^REFERENCE\=(.*)/) {
      our $reference = $1;
    }  
    if ($ln=~/^NGVEC\=(.*)/) {
      our $ngvec = $1;
    }  
	
  }
  close IN;  


  open(IN, "<$lst") or die "Can't open $lst";
  while (<IN>) {
    next if ($_ =~ /^#/);
	$ln = $_;
	chomp($ln);

	if ($ln=~/^SET:((\S+?);(\S+))\s?$/){

	push @sets, $1;
	my @ctrls=split(",", $2);
	my @cases=split(",", $3);
	
	foreach my $ctrl (@ctrls){$samples->{$ctrl}=1;}
	foreach my $case (@cases){$samples->{$case}=1;}
	# only unique keys will be saved
	
    }

  }
  close IN;
  
  
  
}
sub usage {

  my $self = shift;
  my $msg = <<_EOUSAGE_;
  
# CMD:
# perl $0 -ini=ini.txt -lst=list.txt > $0.err 2>&1 
  
# NOTICE:
# 1. the fastq file must be in "??_1.fq" and "??_2.fq" format (decompressed) and placed under ./reads/ dir 
# 2. make sure the ini.txt and list.txt file comply with the below format 
  
# REQUIRED:
ini.txt 
# REFERENCE=<bowtie2 index title>
# NGVEC=<ngvector for ebseq use>
REFERENCE=/export/arrayPRO2/DataBase/mm10/RSEM/reference/mm10.refseq_125polyA
NGVEC=/export/arrayPRO2/DataBase/mm10/RSEM/reference/mm10.refseq_125polyA.transcripts.ngvec

list.txt
# SET:<case1[,case2,...]>;<control1[,control2,...]>
# NOTE: the label should be exactly the readname
SET:minus-minus;minus-plus	
SET:minus-plus;12	
SET:minus-plus;21	
SET:minus-plus;7nM	
SET:minus-plus;CF02	
SET:minus-plus;CN7	
SET:minus-plus;N1	

_EOUSAGE_
;
  return $msg;
}

sub CalculateExpression {
  my $self = shift;
  our $out_CE;
  our $parameters;  
  our $reference;
  our $samples;
  my $item;
  my $r1;
  my $r2;
  my $cmd;
  my $cmd2;
  my @cmds;
  my @cmds2;
  foreach $item (sort keys %{$samples}) {
	if ((-e "$out_CE/$item\.transcript.sorted.bam.bai" ) || (-e "$out_CE/$item\.genome.sorted.bam.bai" )){ print "$item already done mapping!\n"; next;}
    $r1 = "reads/$item\_1.fq";  
    $r2 = "reads/$item\_2.fq";
    # $cmd = "rsem-calculate-expression --bowtie2 -p 20 --paired-end --forward-prob 0 --output-genome-bam --phred33-quals $r1 $r2 $reference $item";   	
	$cmd = "rsem-calculate-expression --bowtie2 -p 20 --paired-end --forward-prob 0  --phred33-quals $r1 $r2 $reference $item";   	
    # $cmd = "rsem-calculate-expression --bowtie2 -p 10 --output-genome-bam --phred64-quals $r1 $reference $item";   
	# remove params (--paired-end --strand-specific --forward-prob 0) to fit BGI's data
	
    # --forward-prob <double>
        # Probability of generating a read from the forward strand of a
        # transcript. Set to 1 for a strand-specific protocol where all
        # (upstream) reads are derived from the forward strand, 0 for a
        # strand-specific protocol where all (upstream) read are derived from
        # the reverse strand, or 0.5 for a non-strand-specific protocol.
        # (Default: 0.5)

    # --strand-specific
        # The RNA-Seq protocol used to generate the reads is strand specific,
        # i.e., all (upstream) reads are derived from the forward strand. This
        # option is equivalent to --forward-prob=1.0. With this option set, if
        # RSEM runs the Bowtie/Bowtie 2 aligner, the '--norc' Bowtie/Bowtie 2
        # option will be used, which disables alignment to the reverse strand
        # of transcripts. (Default: off)		
		
    push (@cmds, $cmd);
    $cmd2 = "mv $item\*  $out_CE";
    push (@cmds2, $cmd2);
  }
  process_cmds_parallel (@cmds);
  process_cmds_parallel (@cmds2);
}

sub GenerateDataMatrix {
  my $self = shift;
  our @controls;
  our @cases;  
  our $out_CE;
  our $out_DM;
  our @sets;
  my ($gene_cmd, $iso_cmd, @cmds);

# # Usage: rsem-generate-data-matrix sampleA.[alleles/genes/isoforms].results sampleB.[alleles/genes/isoforms].results ... > output_name.matrix
# # All result files should have the same file type. The 'expected_count' columns of every result file are extracted to form the data matrix.

  foreach my $set (@sets){
  my @semicolons=split(";", $set); # should be only 2 elements so far
  my $ctrl_str=$semicolons[1];
  my $case_str=$semicolons[0];
  my @ctrls=split(",", $ctrl_str);
  my @cases=split(",", $case_str);  
  my $prefix="${case_str}-vs-${ctrl_str}"; 
  my (@ctrls_gen, @ctrls_iso, @cases_gen, @cases_iso);
  
	foreach my $ctrl_i (@ctrls){
	
	push @ctrls_gen, "${out_CE}/${ctrl_i}.genes.results";
	push @ctrls_iso, "${out_CE}/${ctrl_i}.isoforms.results";
	
	}
	
	foreach my $case_i (@cases){

	push @cases_gen, "${out_CE}/${case_i}.genes.results";
	push @cases_iso, "${out_CE}/${case_i}.isoforms.results";	
	
	}
  
  

	  $gene_cmd = "rsem-generate-data-matrix ".join(" ", @cases_gen)." ".join(" ", @ctrls_gen)."|  sed 's,output/CalculateExpression/,,g' > $out_DM/$prefix\.genMat";
      $iso_cmd  = "rsem-generate-data-matrix ".join(" ", @cases_iso)." ".join(" ", @ctrls_iso)."|  sed 's,output/CalculateExpression/,,g' > $out_DM/$prefix\.isoMat"; 
	  push (@cmds, $gene_cmd);
      push (@cmds, $iso_cmd);	
	  # print "$ctrl_gen\t$case_gen\n";#$ctrl_iso\t$case_iso\n";
	}			
  
  process_cmds_parallel (@cmds);  
  
}

sub ebseq {
  my $self = shift;
  our @controls;
  our @cases;  
  our $out_CE;
  our $out_DM;
  our @sets;  
  our $ngvec;
  our @cmds_2;
  
  # our @ctrldup;
  # our @casedup;
  # our $out_DM;
  our $out_ebseq;
  # our $ngvec;
  my @cmds;
  my $gene_cmd = "";
  my $iso_cmd  = "";
  # my $ctrno;
  # my $caseno;
  # my $idx1;
  # my $idx3;
  my $genctrl_cmd = "";
  my $isoctrl_cmd = "";
  
  
  # for ($idx1 = 0;$idx1 < scalar @cases; $idx1++) {
    # $caseno = $casedup[$idx1];	
	# for ($idx3 = 0;$idx3 < scalar @controls; $idx3++) {	    
	  # $ctrno = $ctrldup[$idx3];  
	  
  foreach my $set (@sets){
  my @semicolons=split(";", $set); # should be only 2 elements so far
  my $ctrl_str=$semicolons[1];
  my $case_str=$semicolons[0];
  my @ctrls=split(",", $ctrl_str);
  my @cases=split(",", $case_str);  
  my $prefix="${case_str}-vs-${ctrl_str}"; 
  my $caseno=scalar(@cases);
  my $ctrno=scalar(@ctrls);  
  
	  $gene_cmd = "rsem-run-ebseq $out_DM/$prefix\.genMat $caseno,$ctrno $out_ebseq/$prefix\.genMat.result";	               
      $iso_cmd  = "rsem-run-ebseq --ngvector $ngvec $out_DM/$prefix\.isoMat $caseno,$ctrno $out_ebseq/$prefix\.isoMat.result"; 
	  $genctrl_cmd = "rsem-control-fdr $out_ebseq/$prefix\.genMat.result 0.05 $out_ebseq/$prefix\.genMat.fdr-0.05.result";
	  $isoctrl_cmd = "rsem-control-fdr $out_ebseq/$prefix\.isoMat.result 0.05 $out_ebseq/$prefix\.isoMat.fdr-0.05.result";
	  push (@cmds, $gene_cmd);
      push (@cmds, $iso_cmd);	
	  push (@cmds_2, $genctrl_cmd);
      push (@cmds_2, $isoctrl_cmd);
	  # print "$gene_cmd\n$genctrl_cmd\n";
	
  }    

  process_cmds_parallel (@cmds);   
  process_cmds_parallel (@cmds_2);   
}


# sub DEGexp {
  # my $self = shift;
  # our @controls;
  # our @cases;
  # our @ctrldup;
  # our @casedup;
  # my $idx1;
  # my $idx3;
  # my $infn;
  # for ($idx1 = 0;$idx1 < scalar @cases; $idx1++) {
	# for ($idx3 = 0;$idx3 < scalar @controls; $idx3++) {	    	  
      # $self->runRscript ($controls[$idx1], $cases[$idx3], $ctrldup[$idx1], $casedup[$idx3]);	  
	# }			
  # }      
# }
# sub runRscript {
  # my $self = shift;
  # my $ctrl = shift;
  # my $case = shift;
  # my $ctrlno = shift;
  # my $caseno = shift;   
  # our $out_DM;  
  # my $genMat = "$out_DM/$case\-$ctrl\.genMat";
  # my $isoMat = "$out_DM/$case\-$ctrl\.isoMat";
  # our $out_DEGseq;
  # my $out_genMat = "$out_DEGseq/$case\-$ctrl\.genMat";
  # my $out_isoMat = "$out_DEGseq/$case\-$ctrl\.isoMat";
  # unless (-e $out_genMat) {system ("mkdir -p $out_genMat");}
  # unless (-e $out_isoMat) {system ("mkdir -p $out_isoMat");}
  
  # my $genrscript = "$out_DEGseq/$case\-$ctrl\.genMat/$case\-$ctrl\.r";
  # my $isorscript = "$out_DEGseq/$case\-$ctrl\.isoMat/$case\-$ctrl\.r";
  
  # my $ctrlvalCol = "";
  # my $casevalCol = "";
  # my $expCol1 = "";
  # my $expCol2 = "";
  # my $count = 1;
  # my $idx = 2;
  # my $cmd;
  # while ($count <= $caseno) {
    # $casevalCol = $casevalCol.$idx.",";
	# $idx = $idx + 1;
    # $count = $count + 1;
  # }
  # chop($casevalCol);  
  # $count = 1;
  # while ($count <= $ctrlno) {
    # $ctrlvalCol = $ctrlvalCol.$idx.",";
	# $idx = $idx + 1;
    # $count = $count + 1;
  # }
  # chop($ctrlvalCol);

  # $count = 1;
  # $idx = 2;
  # while ($count <= $caseno) {
    # $expCol1 = $expCol1.$idx.",";
	# $idx = $idx + 1;
    # $count = $count + 1;
  # }
  # chop($expCol1);  
  
  # $count = 1;
  # $idx = 2;  
  # while ($count <= $ctrlno) {
    # $expCol2 = $expCol2.$idx.",";
	# $idx = $idx + 1;
    # $count = $count + 1;
  # }
  # chop($expCol2);

  # open (OUT, ">$genrscript") or die "Can't create $genrscript\n";
  # print OUT "suppressPackageStartupMessages(library(DEGseq))\n";
  # print OUT "fn <- \"$genMat\"\n";
  # print OUT "case <- readGeneExp(file=fn, geneCol=1, valCol=c($casevalCol))\n";
  # print OUT "ctrl <- readGeneExp(file=fn, geneCol=1, valCol=c($ctrlvalCol))\n";
  # print OUT "layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))\n";
  # print OUT "par(mar=c(2, 2, 2, 2))\n";
  # print OUT "DEGexp(geneExpMatrix1=ctrl, geneCol1=1, expCol1=c($expCol2), groupLabel1=\"$ctrl\", geneExpMatrix2=case, geneCol2=1, expCol2=c($expCol1), groupLabel2=\"$case\", pValue=1e-3, zScore=4, qValue=1e-3, foldChange=2, thresholdKind=4,outputDir=\"$out_genMat\", method=\"MARS\")\n";
  # close OUT;
  # $cmd = "Rscript $genrscript";
  # process_cmd ($cmd);
  # $cmd = "sed 's/value1/$ctrl/' $out_genMat/output_score.txt | sed 's/value2/$case/' > $out_genMat/myoutput_score.txt";
  # process_cmd ($cmd);  
  

  # open (OUT, ">$isorscript") or die "Can't create $isorscript\n";
  # print OUT "suppressPackageStartupMessages(library(DEGseq))\n";
  # print OUT "fn <- \"$genMat\"\n";
  # print OUT "case <- readGeneExp(file=fn, geneCol=1, valCol=c($casevalCol))\n";
  # print OUT "ctrl <- readGeneExp(file=fn, geneCol=1, valCol=c($ctrlvalCol))\n";
  # print OUT "layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))\n";
  # print OUT "par(mar=c(2, 2, 2, 2))\n";
  # print OUT "DEGexp(geneExpMatrix1=ctrl, geneCol1=1, expCol1=c($expCol2), groupLabel1=\"$ctrl\", geneExpMatrix2=case, geneCol2=1, expCol2=c($expCol1), groupLabel2=\"$case\", pValue=1e-3, zScore=4, qValue=1e-3, foldChange=2, thresholdKind=4,outputDir=\"$out_isoMat\", method=\"MARS\")\n";
  # close OUT;
  # $cmd = "Rscript $isorscript";
  # process_cmd ($cmd);
  # $cmd = "sed 's/value1/$ctrl/' $out_isoMat/output_score.txt | sed 's/value2/$case/' > $out_isoMat/myoutput_score.txt";
  # process_cmd ($cmd);  
# }


############################ below are useful subroutines ############################



sub getCurrentDateTime  {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  my $nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
  $year+1900,$mon+1,$mday,$hour,$min,$sec);
  return "[$nice_timestamp]";
}
sub printlog {
  my $msg  = shift;
  open (OUT, ">> my.log") or die "Can't create my.log file \n";
  print OUT "$msg\n";
  close OUT;
}
sub process_cmds_parallel {
  my @cmds = @_;
  my @threads;
  my $msg;
  foreach my $cmd (@cmds) {
  # should only be 2 cmds max
    my $thread = threads->create('process_cmd', $cmd);
    push (@threads, $thread);
  }
                
  my $ret = 0;   
  foreach my $thread (@threads) {
    $thread->join();
    if (my $error = $thread->error()) {
      $msg = "ERROR, thread exited with error $error\n";
      printlog($msg);
      $ret++;
    }
  }
  if ($ret) {
    $msg = "ERROR, $ret threads errored out";
    printlog($msg);
  }
  return;
}
sub process_cmds_serial {
  my @cmds = @_;
  foreach my $cmd (@cmds) {
    process_cmd($cmd);
  }
  return;
}
sub process_cmd {
  my ($cmd) = @_;
  my $start_time = getCurrentDateTime();
  my $ret = system("bash", "-c", $cmd);
  my $end_time = getCurrentDateTime();
  my $msg = "CMD FINISHED: '$cmd' $start_time - $end_time \n";
  if ($ret) {
	chomp $ret;
    $msg = "CMD ERROR: '$cmd' died with ret $ret $start_time - $end_time \n";
    printlog($msg);
  	return;
  }    
  printlog($msg);
  return;
}

1;
