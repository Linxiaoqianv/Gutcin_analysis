mkdir -p ./mapping_result/

# paladin index

paladin index -r3 ./00.data/overall_High.id_0.95_co_0.95.fas

#mapping
for i in SRR*
do
    file=`ls $i | sed -n '1p'`
    base=${file%_clean_fwd.fq.gz}
    cat $i/${base}_clean_fwd.fq.gz $i/${base}_clean_rev.fq.gz > $i/merge.${base}_clean.fq.gz
    paladin align -t 48 ./00.data/overall_High.id_0.95_co_0.95 $i/merge.${base}_clean.fq.gz > ./mapping_result/align_${base}.sam
    samtools view -b ./mapping_result/align_${base}.sam -@ 48 > mapping_result/${base}.bam
    samtools sort ./mapping_result/${base}.bam -o mapping_result/${base}.sort.bam -@ 48
    samtools index ./mapping_result/${base}.sort.bam -@ 48
    samtools idxstats ./mapping_result/${base}.sort.bam > mapping_result/${base}.idxstats.txt
done

#profile
python ./get_count_table.py ./mapping_result/*.idxstats.txt > ./00.data/counts.txt

#analysis
Rscript analysis_code.R