The scripts in this repo were created to analyse the impact of SNVs on splice sites strength. Files:
- 1_SNV_examination - exploring statistics of SNVs from the input VCF files (retrieved from WGS and nascent RNA-seq of Kasumi-1 cell line);
- 2_splicing_sites_selection - creation of set of splice sites from 2 sources: exon-exon junctions previously identified via Kasumi-1 RNA-seq reads mapping and UCSC (intron annotations);
- 3_SNVs_in_splicing_sites - identifying SNVs which are located in splice sites;
- 4_splicing_sites_examination - visualising frequency of bases in each position of splice sites;
- 5_splicing_sites_scoring - assigning scores to reference and altered splice sites and exploring the impact of SNVs in each splice site position on its strength.

PDF versions of each markdown notebook are available.