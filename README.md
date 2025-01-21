# microbial-community-before-and-after-dilution-cultivation
In this study, we diluted raw microbial communities, and analyzed the composition variation in thousands of sub-communities cultured under the same condition.

All the raw sequencing data have been deposited in the NCBI Sequence Read Archive (SRA) database under the BioProject ID PRJNA1066933.The code assumes that the raw data labeled as "TA-raw" and "TA-sub" have already been downloaded to your local computer, and that both VSEARCH and QIIME 2 are installed.

“vsearch.sh” processes the data downloaded under BioProject ID PRJNA1066933 ("TA-raw.fa" and "TA-sub.fa") into a feature table (ASV_table.txt) and representative sequences (ASV.fa).
“qiime.sh” performs analyses such as alpha diversity, beta diversity, and taxonomic annotation using ASV_table.txt and ASV.fa.
“sample-metadata.tsv” contains sample metadata.

Contact: wangjianing@sdu.edu.cn
