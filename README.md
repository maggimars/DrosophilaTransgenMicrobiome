# Drosophila Transgenerational Microbiome Experiment

`source activate qiime2-2017.11`
## import the data ...

  `Csamps=/Users/brisbin/Desktop/carmen2/seqs`
  Check: `echo $Csamps`

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $Csamps \
    --source-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path demux-paired-end.qza

prepare the data for visualization

      qiime demux summarize \
      --i-data demux-paired-end.qza \
      --o-visualization demux.qzv

Choose forward and reverse read cut-off lengths based on Interactive Quality Plots, which are found in the second tab in `demux.qzv`    
  * Forward 280   
  * Reverse 220

## DaDa2 with paired-end reads
`OUTPUTDIR=/Users/brisbin/desktop/carmen2/`

      qiime dada2 denoise-paired \
      --i-demultiplexed-seqs demux-paired-end.qza \
      --output-dir dada2/ \
      --o-representative-sequences rep-seqs \
      --p-trim-left-f 20 \
      --p-trim-left-r 20 \
      --p-trunc-len-f 280 \
      --p-trunc-len-r 220 \
      --p-n-threads 3

check progress `cat /var/folders/s9/xsh5jf0x7y7cgn_sv_7l4b940000gn/T/qiime2-q2cli-err-*`

      qiime feature-table summarize \
      --i-table ./dada2/table.qza \
      --o-visualization ./dada2/table.qzv \
      --m-sample-metadata-file Metadata2.txt

*add taxonomy (Silva 128, 97% 16S)*      

97% taxonomy

Rep_set:
`SILVA97otus=/Users/brisbin/Documents/TonanMaruSequencing/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta`

Taxonomy:
`Tax97=/Users/brisbin/Documents/TonanMaruSequencing/SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_all_levels.txt`


    qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path $SILVA97otus \
    --output-path 97_otus.qza

    qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --source-format HeaderlessTSVTaxonomyFormat \
    --input-path $Tax97 \
    --output-path 97ref-taxonomy.qza

    qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads 97_otus.qza\
    --i-reference-taxonomy 97ref-taxonomy.qza\
    --o-classifier 97classifier.qza

    qiime feature-classifier classify-sklearn \
    --i-classifier 97classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

**filter feature table to include only bacterial sequences**

    qiime taxa filter-table \
    --i-table ./dada2/table.qza \
    --i-taxonomy taxonomy.qza \
    --p-include D_0__Bacteria \
    --o-filtered-table ./dada2/table-bact_only.qza

    qiime feature-table summarize \
    --i-table ./dada2/table-bact_only.qza \
    --o-visualization ./dada2/table-bact_only.qzv \
    --m-sample-metadata-file Metadata2.txt

**Create Bacteria-Only Taxonomy Barplots visualization for all samples**      

    qiime taxa barplot \
    --i-table ./dada2/table-bact_only.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file Metadata2.txt \
    --o-visualization taxa-bar-plots-bact_only.qzv

remove media samples:

    qiime feature-table filter-samples \
    --i-table ./dada2/table-bact_only.qza \
    --m-metadata-file Metadata2.txt \
    --p-where "SampleType='Gut'" \
    --o-filtered-table ./dada2/gut-table-bact_only.qza

    qiime feature-table summarize \
    --i-table ./dada2/gut-table-bact_only.qza \
    --o-visualization ./dada2/gut-table-bact_only.qzv \
    --m-sample-metadata-file Metadata2.txt

**Make a phylogenetic tree for phylogenetic diversity analysis (Unifrac and Weighted Unifrac)**   

    qiime alignment mafft \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza

    qiime alignment mask \
    --i-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza

    qiime phylogeny fasttree \
    --i-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza

    qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree

**phylogenetic diversity metrics**

    qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --m-metadata-file Metadata2.txt \
    --i-table ./dada2/gut-table-bact_only.qza \
    --p-sampling-depth 1000 \
    --output-dir core-metrics-phylogenetic


**Filter distance matrices**    
FO     

      qiime diversity filter-distance-matrix \
      --i-distance-matrix core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza \
      --m-metadata-file Metadata2.txt \
      --p-where "Generation='F0'" \
      --o-filtered-distance-matrix core-metrics-phylogenetic/F0-weighted_unifrac_distance_matrix.qza   
F1

    qiime diversity filter-distance-matrix \
    --i-distance-matrix core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file Metadata2.txt \
    --p-where "Generation='F1'" \
    --o-filtered-distance-matrix core-metrics-phylogenetic/F1-weighted_unifrac_distance_matrix.qza    

F2

      qiime diversity filter-distance-matrix \
      --i-distance-matrix core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza \
      --m-metadata-file Metadata2.txt \
      --p-where "Generation='F2'" \
      --o-filtered-distance-matrix core-metrics-phylogenetic/F2-weighted_unifrac_distance_matrix.qza    

F3

      qiime diversity filter-distance-matrix \
      --i-distance-matrix core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza \
      --m-metadata-file Metadata2.txt \
      --p-where "Generation='F3'" \
      --o-filtered-distance-matrix core-metrics-phylogenetic/F3-weighted_unifrac_distance_matrix.qza            

**PERMANOVA pair-wise by Treatment within generations**

F0

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-phylogenetic/F0-weighted_unifrac_distance_matrix.qza \
    --m-metadata-file Metadata2.txt \
    --m-metadata-category Treatment \
    --o-visualization core-metrics-phylogenetic/F0-weighted_unifrac-Treatment-significance.qzv \
    --p-pairwise

F1

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-phylogenetic/F1-weighted_unifrac_distance_matrix.qza \
    --m-metadata-file Metadata2.txt \
    --m-metadata-category Treatment \
    --o-visualization core-metrics-phylogenetic/F1-weighted_unifrac-Treatment-significance.qzv \
    --p-pairwise

F2

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-phylogenetic/F2-weighted_unifrac_distance_matrix.qza \
    --m-metadata-file Metadata2.txt \
    --m-metadata-category Treatment \
    --o-visualization core-metrics-phylogenetic/F2-weighted_unifrac-Treatment-significance.qzv \
    --p-pairwise

F3

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-phylogenetic/F3-weighted_unifrac_distance_matrix.qza \
    --m-metadata-file Metadata2.txt \
    --m-metadata-category Treatment \
    --o-visualization core-metrics-phylogenetic/F3-weighted_unifrac-Treatment-significance.qzv \
    --p-pairwise

## export data for phyloseq

**Taxonomy**        

    qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

Drag and drop `taxonomy.qzv` in view.qiime2.org window. In the top left of the window there are 'export options' - click export to TSV. The file will automatically download as `metadata.csv`. Move the file to Carmen2 directory and rename `Taxonomy.csv`

**Phylogenetic Tree**    

    qiime tools export \
    ./rooted-tree.qza \
    --output-dir ./

outPut = `tree.nwk`


**OTU table**     
    qiime tools extract \
    ./dada2/table.qza \
    --output-dir ./

convert `.biom` file to `.txt` file for R

    biom convert -i ./ExtractedFeatureTable/data/feature-table.biom -o ./feature-table.txt --to-tsv

Open text file and format - remove `#Constructed from biom file` line and the `#` in front of the first column name
