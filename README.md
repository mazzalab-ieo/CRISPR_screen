# CRISPR_screen

CRISPR_scripts:

**1.mageck_count.sh**  —>  This script is the first one of the our pipeline and it used to run MAGeCK count for processing CRISPR screening data. It takes as input FASTQ files and generate as output a matrix collecting all read counts saved in the specified directory.  It specifies the number of bases to trim from the 5’ end of each read in order to remove adapter sequences or low-quality bases and also specifies the sgRNA list used on the CRISPR screen. The output is a matrix collecting all read counts. 

**2.crispr_gw.R** —> This is the second script of our pipeline and it performs counts normalisation, computation of fold changes, correction for gene-independent effects and results visualization with “CRISPRcleanR” package.

**3.mageck_mle.sh** —> This it the third script of our pipeline, it takes the normalised counts from the previous script and runs MAGeCK MLE which performs a maximum likelihood analysis of gene essentialities. It needs a design matrix file which contains information about the experimental setup.

**4.mageck_flute.Rmd** —> This is the forth and final script of our pipeline which runs MAGeCKFlute for accurate identification of gene hits and associated pathways. It takes as input the count summary file obtained from the mageck_count.sh to perform quality control and the gene summary file from mageck_mle.sh for the downstream analysis. 
