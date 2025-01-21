#Use the obtained “ASV_table.txt” and “ASV.fa” files from the VSEARCH for further analysis with QIIME2. sample-metadata.tsv contains sample metadata.

#Data preprocessing.
mkdir qiime2
cd qiime2
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}' ASV.fa >ASV2.fa

#Import the data into QIIME 2.
conda activate qiime2-2020.2
biom convert -i ASV_table.txt -o ASV_table.biom --table-type="OTU table" --to-json
qiime tools import --input-path ASV_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path XM11-YJ12-table.qza
qiime tools import --input-path ASV2.fa --output-path ASV2.qza --type 'FeatureData[Sequence]'

#View the statistics of the feature table.
qiime feature-table summarize --i-table XM11-YJ12-table.qza --o-visualization XM11-YJ12-table.qzv --m-sample-metadata-file sample-metadata.tsv

#Filter samples with the reads at least 5000.
qiime feature-table filter-samples --i-table XM11-YJ12-table.qza --p-min-frequency 5000 --o-filtered-table XM11-YJ12-5000-filtered-table.qza

#Count the representative sequences
qiime feature-table tabulate-seqs --i-data ASV2.qza --o-visualization ASV2-rep-seqs.qzv

#tree and diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ASV2.qza --o-alignment ASV2-aligned-rep-seqs.qza --o-masked-alignment ASV2-masked-aligned-rep-seqs.qza --o-tree ASV2-unrooted-tree.qza --o-rooted-tree ASV2-rooted-tree.qza
qiime diversity core-metrics-phylogenetic --i-phylogeny ASV2-rooted-tree.qza --i-table XM11-YJ12-table.qza --p-sampling-depth 4900 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/observed_otus_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/observed-otus-group-significance.qzv

#Conduct taxonomic analyses.
qiime feature-classifier classify-sklearn --i-classifier /mnt/sdb/wjn/96kong/training-feature-classifiers/classifier.qza --i-reads ASV2.qza --o-classification ASV2-taxonomy.qza
qiime metadata tabulate --m-input-file ASV2-taxonomy.qza --o-visualization ASV2-taxonomy.qzv
qiime taxa barplot --i-table XM11-YJ12-table.qza --i-taxonomy ASV2-taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization XM11-YJ12-taxa-bar-plots.qzv











