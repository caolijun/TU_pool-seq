#! /usr/bin/perl -w
use strict;

system('mkdir -p /work/home/caolijun/tu_pool_seq/1rawdata/clean/scripts');
my($fastp,$core,$q,$suffix,$IN,$OUT);
$IN="/work/home/caolijun/tu_pool_seq/1rawdata/tu_2017";
$suffix=".R1.fastq.gz";
$OUT='/work/home/caolijun/tu_pool_seq/1rawdata/clean';
$fastp='/work/share/kdy_caolijun/software/fastp';
$core =6;
$q="xhacnormalb";


my @list1s=`find $IN -name "*$suffix"`;
my(@lists);
push @lists,@list1s;
print "$#lists\n";

&createJob(@lists);
&submitJob(@lists);


sub submitJob{
        my @arr=@_;
        foreach (@arr){
           chomp;
           my @arrays = split/\//;
           my $sample = $arrays[-1];
           print $sample,"\n";
           $sample =~s/$suffix//; 
           my $Jobs=`squeue | grep "fastp" | wc -l`;
           while( $Jobs > 98 ){
              print "JOB Remain $Jobs\n";
              sleep 98;
              $Jobs=`squeue | grep "fastp" | wc -l`;
           }
           system("sbatch $OUT/scripts/$sample.pbs >$OUT/scripts/$sample.log");
        }
}
sub createJob{
        my @arr=@_;
        foreach (@arr){
                chomp;
                my $R1 = $_;
                my @arrays = split/\//,$R1;
                my $sample = $arrays[-1];
                $sample =~s/$suffix//; 
                my $R2=$_;
                $R2=~s/.R1/.R2/;
                open OUT,">$OUT/scripts/$sample.pbs";
                print OUT <<SET;
#!/bin/bash
#SBATCH -J fastp_$sample
#SBATCH -p $q
#SBATCH -n $core
#SBATCH -N 1
#SBATCH --mem=18g
#SBATCH -o out$sample.%j
#SBATCH -e err$sample.%j
date
#callsnp

$fastp -w $core -i $R1 -I $R2 -o $OUT/$sample\_R1.fastq.gz -O $OUT/$sample\_R2.fastq.gz
rm $R1 $R2
date
SET
        close OUT;
        }
}
