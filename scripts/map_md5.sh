while read sample read1 read2 
do
    #paste md5sums/${sample}_md5sum_r1.tsv | awk -F'[\t/]' '{print  $NF"\tfastq\t"$1}'  
    #paste md5sums/${sample}_md5sum_r2.tsv | awk -F'[\t/]' '{print  $NF"\tfastq\t"$1}'  
    f1=`echo $read1 | awk -F'/' '{print $NF}'`
    f2=`echo $read2 | awk -F'/' '{print $NF}'`
    echo -e "$f1\t$f2"
    
    
done < metadata/md5_gen.tsv
