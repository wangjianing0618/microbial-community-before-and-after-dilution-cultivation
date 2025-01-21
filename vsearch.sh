mkdir -p temp result
#Place the downloaded files (TA-sub.fa and TA-raw.sub) into the temp folder. “TA-sub.fa” contains 16S rRNA gene sequencing data from the communities after dilution-cultivation, while “TA-raw.fa” represents the original communities.
cat temp/TA-sub.fa temp/TA-raw.fa >> temp/qc.fa

vsearch --fastx_filter temp/qc.fa \
    	--fastq_stripleft 19 --fastq_stripright 18 \
    	--fastaout temp/filtered.fa

#Calculate frequency of non-redundancy reads
    vsearch \
        --derep_fulllength temp/filtered.fa \
    	--relabel Uni --minuniquesize 8 --sizeout \
    	--output temp/uniques.fa 

#Denoise by unoise3
    usearch -unoise3 temp/uniques.fa \
        -zotus temp/Zotus.fa

#Rename to ASV
    awk 'BEGIN {n=1}; />/ {print ">ASV_" n; n++} !/>/ {print}' temp/Zotus.fa \
        > result/ASV.fa

	vsearch --usearch_global temp/filtered.fa \
	    --db result/ASV.fa \
        --otutabout temp/ASV_table.txt \
        --id 0.97
#Use the obtained “ASV_table.txt” and “ASV.fa” files for further analysis with QIIME2.

    
