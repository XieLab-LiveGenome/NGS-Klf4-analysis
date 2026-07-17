echo \#!/bin/bash
echo
echo runfolder=\"/data2/xiongl2/Projects/Rawdata/$3/Samfiles_spikein/$2\"
echo mkdir \$runfolder
echo cd \$runfolder
echo
echo "/data2/xiongl2/Softwares/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --genomeDir /data2/xiongl2/Projects/annotation/spikein --readFilesIn \$(ls /data2/xiongl2/Projects/Rawdata/$3/Fastqfiles/$2_*.fastq.gz) --runThreadN $1 --outFilterMismatchNmax 2 --outFilterMultimapScoreRange 0 --alignIntronMax 500000 --readFilesCommand zcat &>output-$2-run.STAR"
echo
echo "samtools view -b -h -q 7 -f 2 -@ $1 Aligned.out.sam |samtools sort -@ $1 -o $2.bam"
echo "samtools index $2.bam"
echo
echo rm -f Aligned.out.sam
echo