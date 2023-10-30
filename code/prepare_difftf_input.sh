input=/home/PERSONALE/nicola.balboni2/agc1_atac/data/006_difftf/input
index=${input}/referenceGenome/bt_index/mm10
cd $input/new_fastq

# Align reads
for file in *_R1_001.fastq.gz
do
root=`basename $file _R1_001.fastq.gz`
out=${input}/bowtie_bam/out_${root}.txt
echo "Aligning sample $root"
file1=${root}_R1_001.fastq.gz
file2=${root}_R2_001.fastq.gz
bowtie2 -x $index -1 $file1 -2 $file2 -p 24 --met-file $out | samtools view -@ 24 -b -S -h - > ${input}/bowtie_bam/${root}.bam

echo "Sorting sample $root"
samtools sort -m 5G -@ 24 -O BAM -o ${input}/bowtie_bam/${root}.sorted.bam ${input}/bowtie_bam/${root}.bam

echo "Calling peaks on sample $root"
macs2 callpeak -t ${input}/bowtie_bam/${root}.sorted.bam -f BAMPE \
--outdir ${input}/peaks/ --name $root \
-g mm -q 0.05
done

# Align rnaseq reads
# cd $input/rnaseq_fastq
# for file in *_R1_001.fastq.gz
# do
# root=`basename $file _S1_L001_R1_001.fastq.gz`
# out=${input}/rnaseq_bam
# echo "Aligning sample $root"
# file1=${root}_S1_L001_R1_001.fastq.gz
# file2=${root}_S1_L001_R2_001.fastq.gz
# bowtie2 -x $index -1 $file1 -2 $file2 -p 18 --met-file $out/out_${root}.txt | samtools view -@ 18 -b -S -h - > ${input}/rnaseq_bam/${root}.bam
# 
# echo "Sorting sample $root"
# samtools sort -m 5G -@ 18 -O BAM -o ${out}/${root}.sorted.bam ${out}/${root}.bam
# done
# 
# echo "counting reads in sample $root"
# featureCounts -T 18 -s 1 -g gene_id -p -a ${input}/referenceGenome/mm10.gtf.gz -o ${out}/agc1_mm10.counts ${out}/*.sorted.bam


