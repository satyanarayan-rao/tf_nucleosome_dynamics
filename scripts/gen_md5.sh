while read sample read1 read2 
do
    echo "md5sum $read1 > md5sums/${sample}_md5sum_r1.tsv" 
    echo "md5sum $read2 > md5sums/${sample}_md5sum_r2.tsv" 
done < metadata/md5_gen.tsv
