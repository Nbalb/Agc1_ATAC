cd $agc1v2

for rname in data/new_bam/*.sorted.bam

do

outname=`basename $rname .sorted.bam`
MACS2 callpeak -t $rname -f BAMPE \
--outdir $agc1/data/peaks/ --name $outname \
-g mm -q 0.05

done
