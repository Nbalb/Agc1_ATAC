# Run from the project folder 
#!/bin/sh
for rname in data/001/bam/*.sorted.bam

do

outname=`basename $rname .sorted.bam`
echo "Running macs2 callpeak on sample $outname"
macs2 callpeak -t $rname -f BAMPE \
--outdir data/001/peaks/ --name $outname \
-g mm -q 0.05

done
