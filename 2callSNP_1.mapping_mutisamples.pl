########################################################################
#	File Name: globBWAMEM.pl
#	> Author: QiangGao
#	> Mail: qgao@genetics.ac.cn 
#	Created Time: Fri 11 Nov 2016 11:22:41 AM CST
#	Modified time: Thurs 9 Nov 2023 13:03 PM CST
#	Modified Author: LijunCao
#########################################################################

############nohup to submit the mission


#!/usr/bin/perl
######data setting
use strict;
my $script_dir="/work/home/caolijun/tu_pool_seq/2popoolation/scripts";###including getNM.pl pick_up_mapping_quality.pl
my $REF="/work/home/caolijun/tu_pool_seq/2popoolation/genome/tu.genome.ipm_v2.fasta";
my $input_raw="/work/home/caolijun/tu_pool_seq/1rawdata/clean/clean2"; ###including all samples in the path, and withou sub-folds
my $output_data="/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data"; 
my $output_log="/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_log"; 
my $prefix="tu";
my $suffix="_R1.fastq.gz"; ##change to your data's suffix
my $queue="xhacnormalb";
my $CORE=40;
my $pe=150;
&help unless defined ($prefix && $REF && $input_raw && $suffix);
#######default software diray
my$samtool_dir="/work/share/kdy_caolijun/software/samtools/samtools-1.14/";
my$bwa2='/work/share/kdy_caolijun/software/bwa-mem2-2.2.1/bwa-mem2';
my$picard_dir="/work/share/kdy_caolijun/software/";
my$seqkit_dir="/work/share/kdy_caolijun/software/";
my$Qualimap_dir="/work/share/kdy_caolijun/software/qualimap_v2.2.1";
if ($ARGV[0] ne "yes" && $ARGV[0] ne "YES"){ &help;die;}
sub help{
print "\n";
print "Mission Prefix		= $prefix\n";
print "suffix			= $suffix\n";
print "queue 			= $queue\n";
print "CPU   			= $CORE\n";
print "length			= PE150\n";
print "REF  			= $REF\n";
print "input fastq dir		= $input_raw\n";
print "script			= $script_dir\n";
print "samtools 		= $samtool_dir\n";
print "bwa			= $bwa2\n";
print "picard			= $picard_dir\n";
print "seqkit			= $seqkit_dir\n";
print "qualimap		= $Qualimap_dir\n\n";

print "if all parameter is OK, please type \n\n";
print "\t nohup perl $0 yes & \n\n";
}

#################Start
my @file=`find $input_raw -name "*$suffix" `;


if(-e $output_data){
	
}else{
	my $cc=`mkdir -p $output_data`;
}

if(-e $output_log){
	
}else{
	my $cc=`mkdir -p $output_log`;
}
my$bwa_mem=$CORE*3;
$bwa_mem.="g";
my$name;
my$total_number=scalar(@file);
open(USED,">used.data");
foreach(@file){
	chomp $_;
	($name)=$_=~/.*\/(.*?)$suffix/;
        print $name;
	my $R2=$_;
	$R2=~s/_R1/_R2/;                   #######need change 
	open(OUT,">$output_log/$name.map.slurm");
print OUT <<EOF;
#!/bin/bash
#SBATCH -J $prefix\_bwa$name
#SBATCH -p $queue
#SBATCH -n $CORE
#SBATCH -N 1
#SBATCH --mem=$bwa_mem
#SBATCH -o $output_log/$name.slurm.log
#SBATCH -e $output_log/$name.slurm.err

export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.131-11.b12.el7.x86_64/
export CLASSPATH=.:\$JAVA_HOME/lib:\$JAVA_HOME/jre/lib
export PATH=\$JAVA_HOME/bin:\$JAVA_HOME/jre/bin:\$PATH

### Map data using bwa-mem2
$bwa2 mem -t $CORE -R \"\@RG\\tID:$name\\tLB:$name\\tSM:$name\\tPL:illumina\\tPU:$name\" $REF $_ $R2| $samtool_dir/samtools view -@ $CORE -bS -o $output_data/$name.source.bam
$samtool_dir/samtools flagstat $output_data/$name.source.bam  >$output_data/$name.source.mapinfo
$samtool_dir/samtools sort -m 3g -@ $CORE $output_data/$name.source.bam -o $output_data/$name.noq30.sorted.bam
rm $output_data/$name.source.bam
$samtool_dir/samtools stats $output_data/$name.noq30.sorted.bam > $output_data/$name.noq30.sorted.bam.stat
perl $script_dir/getNM.pl $output_data/$name.noq30.sorted.bam > $output_data/$name.noq30.sorted.bam.NM
$samtool_dir/samtools view -@ $CORE -h -b -q30 $output_data/$name.noq30.sorted.bam > $output_data/$name.sorted.bam
rm $output_data/$name.noq30.sorted.bam
java -Xmx120G -jar $picard_dir/picard.jar MarkDuplicates I=$output_data/$name.sorted.bam O=$output_data/$name.sorted.rmdup.bam CREATE_INDEX=true REMOVE_DUPLICATES=true M=$output_data/$name.marked_dup_metrics.txt
###########
######
rm $output_data/$name.sorted.bam
EOF
  close OUT;
  next if(-e "$output_data/$name.source.mapinfo");
  my $cc=`sbatch < $output_log/$name.map.slurm >> bwa.monitor.list`;
  print "run $name\n";
}
close USED;
#####Genome length
my$ref_length=`$seqkit_dir/seqkit stat ${REF} |tail -1 |awk '{OPS="\s+";print \$5}'`;
$ref_length =~ s/,//g;
my$now_path=`pwd`;

sleep 60;
my$mapinfo_file=`ls $output_data/*mapinfo |wc -l `;
chomp($mapinfo_file);
until($mapinfo_file == $total_number){
        sleep 300;
        print "circle is on .\n";
        $mapinfo_file=`ls $output_data/*mapinfo |wc -l`;
        chomp ($mapinfo_file);
        print "$total_number\t $mapinfo_file\t .\n";
}
my$mapinfo_all=`find $output_data -name "*mapinfo"`;
my@mapinfo_single=split(/\n/,$mapinfo_all);
my$tem_count=0;
for(my$e=0;$e<scalar(@mapinfo_single);$e++){
        if (-z $mapinfo_single[$e]){
        $tem_count++;
        }
}
until($tem_count ==0){
        sleep 60;
        $tem_count=0;
        for(my$e=0;$e<scalar(@mapinfo_single);$e++){
        if (-z $mapinfo_single[$e]){
        $tem_count++;
                }
        }
}

if ($tem_count != 0 ){die "$tem_count\tmapinfo is zero file,error.\n";}
`rm $output_data/*source.bam`;
#####check sorted.bam be killed or not
my$kill=`grep -E "kill|cound't|out-of-memory" $output_log/*err `;
my$kill_mission=`grep -E "kill|cound't|out-of-memory" $output_log/*err |wc -l`;
if ($kill_mission != 0 ){
print "$kill\n";
die "some bwa mission was killed, there are $kill_mission mission was killer.\n";}

######pick up
`perl $script_dir/pick_up_mapping_quality.pl -input $output_data -outdir $output_data -pe $pe -size $ref_length`;
print "pick up is complete.\n";
