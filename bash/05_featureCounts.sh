featureCounts \
-p -T 4 -s 2 -t exon -g gene_id \
-a /mnt/d/rnaseq_data/reference/Oryza_sativa.IRGSP-1.0.63.gtf \
-o /mnt/d/rnaseq_data/counts/gene_counts.txt \
/mnt/d/rnaseq_data/aligned/*_sorted.bam
