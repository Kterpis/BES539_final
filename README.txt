This script is currently set up to run on University of Rhode Island's cluster, Bluewaves. 
The path at line 116 is hardcoded. You will need to change it to your own working directory path. 
You will also need to download the adapters file and set the path to the adapters at lines 46, 58, 70, and 82.

This pipeline starts with downloading the data from the SRA, and renaming the data. 
The transcripts are run through FastQC, trimmed, and assembled.
The headers of each file are manipulated so it is easier to see which header is from which species. 
The transcripts are finally run through OrthoFinder to find single copy orthologs.
After this script has run to completion, download the file called "Statistics_PerSpecies.tsv" located in the Comparative_Genomics_Statistics subdirectory of the Results[date]_1 directory generated as part of the orthofinder output
This file uses the supplied R script to visualize the results
Line 4 of the R script is hardcoded to the location the file was downloaded to. 