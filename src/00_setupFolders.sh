#!/bin/bash
# This script will create the necesarry folders and download the corresponding data from ArrayExpress

mkdir -p ../data/counts/
cd ../data/counts
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8932/E-MTAB-8932.processed.1.zip
unzip E-MTAB-8932.processed.1.zip
