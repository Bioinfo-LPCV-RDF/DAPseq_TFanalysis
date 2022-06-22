
# ---- setting work environment (to be change depending on your (super)computer system)-----
# -- request node on the super computer 
# salloc -c 6 -p big --constraint=sl7 --mem=15g --chdir=$PWD
# -- log to the node 
# srun --jobid=$SLURM_JOB_ID --pty /bin/bash
# -- activate the conda environment
# conda activate DAPseqEASY

main_dir=/home/312.6-Flo_Re/312.6.1-Commun/Romain/book_DAPseqProtocol/conda_bioinfo 
#main_dir=/path/to/your/working/directory #replace with path to your working directory
meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip # v4.12.0 # change to your MEME path
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme # v4.12.0

mspc_mk=/home/312.6-Flo_Re/312.6.1-Commun/lib/MSPC/mspc/mspc # change to your MSPC path

# prepare out directories for each steps
data_dir=$main_dir/data
clean_data_dir=$main_dir/clean_data
mapping_dir=$main_dir/mapping
peaks_dir=$main_dir/peaks
motif_dir=$main_dir/motif
roc_dir=$main_dir/roc
mkdir -p $data_dir $clean_data_dir $mapping_dir $peaks_dir $motif_dir $roc_dir

# general  settings
seed=7539293
threads=6

# format and index the reference genome (here TAIR10 Arabidopsis thaliana) 
if [ ! -f $data_dir/Athaliana_index/at.1.bt2 ]
then
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz 
mv TAIR10_chr_all.fas.gz $data_dir
gunzip $data_dir/TAIR10_chr_all.fas.gz
printf "ChrC\nChrM\n" > ids.txt
awk 'BEGIN{while((getline<"ids.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $data_dir/TAIR10_chr_all.fas > $data_dir/tair10.fas
mkdir $data_dir/Athaliana_index
bowtie2-build --threads $threads $data_dir/tair10.fas $data_dir/Athaliana_index/at
rm ids.txt
fi
PATH_TO_BOWTIE_INDEX="$data_dir/Athaliana_index/at"
Genome_length=120000000 # for peak calling (you need to change this value to fit with your reference genome assembly)
genome="$data_dir/tair10.fas"

# - - - - - - - - - - - - Download toy data - - - - - - - - -
# bHLH BIM2 (2 rep colamp)
control=SRR2926069
TFname="BIM2"
sample1=SRR2926982
sample2=SRR2926983


echo "download fastq data"
do_download () {
prefetch $1 --max-size 30g -O $data_dir
#fasterq-dump $1 -O $data_dir -t $TMPDIR --split-3 -p -e $threads
fastq-dump $1 -O $data_dir --split-3
gzip -c $data_dir/${1}_1.fastq > $data_dir/${1}_1.fastq.gz
gzip -c $data_dir/${1}_2.fastq > $data_dir/${1}_2.fastq.gz
}
#do_download $control
#do_download $sample1
#do_download $sample2
echo "donwload has finished"


echo "quality check and cleaning"
do_QC_cleaning () {
# run FAstQC to check reads quality
fastqc -t $threads -o $data_dir $data_dir/${1}_1.fastq.gz $data_dir/${1}_2.fastq.gz
# run NGmerge in adapter-mode to remove adapters
NGmerge -a -1 $data_dir/${1}_1.fastq.gz -2 $data_dir/${1}_2.fastq.gz -q 34 -o $clean_data_dir/${1}_adapterRm -n $threads
#fastqc -t $threads -o $clean_data_dir $clean_data_dir/${1}_adapterRm_1.fastq.gz $clean_data_dir/${1}_adapterRm_2.fastq.gz
# run Trimmomatics to quality trim reads
trimmomatic PE $clean_data_dir/${1}_adapterRm_1.fastq.gz $clean_data_dir/${1}_adapterRm_2.fastq.gz $clean_data_dir/${1}_final_forward_paired.fastq.gz $clean_data_dir/${1}_final_forward_unpaired.fastq.gz $clean_data_dir/${1}_final_reverse_paired.fastq.gz $clean_data_dir/${1}_final_reverse_unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# run FastQC on the trimmed reads
fastqc -t $threads -o $clean_data_dir $clean_data_dir/${1}_final_forward_paired.fastq.gz $clean_data_dir/${1}_final_forward_unpaired.fastq.gz $clean_data_dir/${1}_final_reverse_paired.fastq.gz $clean_data_dir/${1}_final_reverse_unpaired.fastq.gz
}
#do_QC_cleaning $control
#do_QC_cleaning $sample1
#do_QC_cleaning $sample2
echo "quality check and cleaning has finished"


echo "reads mapping with bowtie2"
# mapping function
do_mapping () {
	bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 $clean_data_dir $clean_data_dir/${1}_final_forward_paired.fastq.gz -2 	$clean_data_dir $clean_data_dir/${1}_final_reverse_paired.fastq.gz -S $mapping_dir/${1}.sam --dovetail -p $threads 2>$mapping_dir/${1}_log.txt
	
	samtools view -Sh $mapping_dir/${1}.sam | \
		grep -e "^@" -e 'XM:i:[012][^0-9]' | awk '$1~/@/ || $5>30 {print $0}' | grep -v "XS:i:" > $mapping_dir/${1}.filtered.sam

	echo "Step 3: SAM to BAM conversion"
	samtools view -Sh -b $mapping_dir/${1}.filtered.sam > $mapping_dir/${1}.filtered.bam
	samtools sort -o $mapping_dir/${1}.filtered.sorted.bam $mapping_dir/${1}.filtered.bam
	echo "Step 4: remove PCR duplicates";
	samtools rmdup -s $mapping_dir/${1}.filtered.sorted.bam $mapping_dir/${1}.filtered.sorted.nodup.bam
	echo "Step 5: BAM indexing";
	samtools index $mapping_dir/${1}.filtered.sorted.nodup.bam
	
}
# launch mapping for the DAP replicates and the control input DNA
do_mapping $sample1
do_mapping $sample2
do_mapping $control
echo "reads mapping has finished"

echo "peaks calling"

if [ ! -f ${peaks_dir}/control.bed ]
then
#  BAM to BED formating for the control mapping results
bamToBed -i $mapping_dir/${control}.filtered.bam > $peaks_dir/control.bed
sortBed -i $peaks_dir/control.bed > $peaks_dir/control.bed.sorted
mv $peaks_dir/control.bed.sorted $peaks_dir/control.bed
fi

# function for peak calling per replicate
do_macs2 () {
#convert to BED
bamToBed -i $mapping_dir/${1}.filtered.bam > $peaks_dir/${1}.bed
sortBed -i $peaks_dir/${1}.bed > $peaks_dir/${1}.bed.sorted
mv $peaks_dir/${1}.bed.sorted $peaks_dir/${1}.bed
#PATH_TO_MACS=$(dirname $(readlink -f $(which macs2)))
log=$peaks_dir/${1}_log.txt
echo "Started MACS..."
macs2 callpeak -t $peaks_dir/${1}.bed -c $peaks_dir/control.bed -B -f BEDPE -n $peaks_dir/${1} -g $Genome_length --call-summits -q 0.00001 >> $log 2>&1;
echo "MACS2 has finished"
#fragment_length=$(cat $log | grep "fragment size = " | awk -v OFS="\t" '{print $12}')
}
do_macs2 $sample1
do_macs2 $sample2

#MSPC
echo "start MSPC"

list_peaks_mspc=("$peaks_dir/${sample1}_peaks.narrowPeak" "$peaks_dir/${sample2}_peaks.narrowPeak")
if [ ! -f $peaks_dir/$TFname/ConsensusPeaks.bed ]
 then
 echo "${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c 100% -o $peaks_dir/$TFname -d 1"
 ${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c 100% -o $peaks_dir/$TFname -d 1
fi
echo "MSPC has finished"
echo "peaks calling has finished"

# - - - - motif analysis and PWM reconstruction
calc(){
awk "BEGIN { print "$*" }"; 
}

 
all_peaks=$peaks_dir/$TFname/ConsensusPeaks_mspc_peaks.txt
learning_size=600
cat  $all_peaks| sort -k6,6nr | head -n $learning_size | cut -f -3 > $motif_dir/training_peaks.bed # rank peaks according to Chisqdv and keep th best 600

bedtools getfasta -fi $genome -bed $motif_dir/training_peaks.bed -fo $motif_dir/training_peaks.fas

echo "meme"

#-meme-pal
#meme-chip -oc $motif_dir/meme -minw 8 -maxw 10 -seed $seed $motif_dir/training_peaks.fas

$meme_prog -oc $motif_dir/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*1000) -meme-minw 8 -meme-maxw 10 -meme-nmotifs 1 -dreme-m 0 -noecho $motif_dir/training_peaks.fas -seed $seed

$meme2meme $motif_dir/meme/meme_out/meme.txt > $motif_dir/meme/meme_out/meme_mini.txt
path_to_meme_mini=$motif_dir/meme/meme_out/meme_mini.txt
bash meme2pfm.sh $path_to_meme_mini $TFname 
cp $motif_dir/meme/meme_out/Motif_MEME_seperateFiles/Motif_1.pfm $motif_dir/${TFname}.pfm

# ROC analysis
#test set
cat  $all_peaks| sort -k6,6nr | tail -n +$learning_size | head -n 5000 | cut -f -3 > $roc_dir/test_peaks.bed
bedtools getfasta -fi $genome -bed $roc_dir/test_peaks.bed -fo $roc_dir/test_peaks.fas

meme_shuffle=/home/prog/meme/meme_4.12.0/bin/fasta-dinucleotide-shuffle
$meme_shuffle -f $roc_dir/test_peaks.fas -s $seed -c 1 > $roc_dir/test_peaks_SHUFFLED.fas

matrice=$motif_dir/${TFname}.pfm
mkdir -p $roc_dir/scores
python scores.py -m ${matrice} -f $roc_dir/test_peaks.fas -o ${roc_dir}/scores/
python scores.py -m ${matrice} -f $roc_dir/test_peaks_SHUFFLED.fas -o ${roc_dir}/scores/
paste <(sort -u -k1,1 -k8,8nr $roc_dir/scores/test_peaks.fas.scores | sort -u -k1,1 | awk '{print $8}'  | sort -nr )  <(sort -u -k1,1 -k8,8nr $roc_dir/scores/test_peaks_SHUFFLED.fas.scores | sort -u -k1,1 | awk '{print $8}' | sort -nr ) | awk '{print $0}' >"$roc_dir/scores/tab_pfm.tsv"

echo "plot ROC"
Rscript ROC.R $roc_dir

echo "finito"
















echo "finito"