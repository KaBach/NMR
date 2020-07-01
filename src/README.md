# Scripts to reproduce the analysis

The easiest way to reproduce the analysis is to use the bash script to download the data and reproduce the figures with NMR.R and NMR_And_MM.Rmd
- 00_setupFolders.sh is a convenience script to download the processed data from array express. 
- NMR.R contains the scripts to compute differential expression between various comparisons.
- NMR_And_MM.Rmd contains the code to produce the PCA plot in the paper.

## Additional scripts
These don't need need to be rerun and are just here for completeness
- GeneIDMappingCleanUp.R contains code to clean up the entrez gene id mapping
- 01_prepareNMR.R processes the output from STAR
