######################################################################################
# CRISPR/Cas9 HT NGS pipeline - MAGeCK v.0.5.9
# Author: Gianluca Vozza, bioinformagician
# Version: 1.0
###################################################################################### 
# D. Trastulli experiment                                                            #
#                                                                                    #
#                                                                                    #
######################################################################################
source activate /hpcnfs/scratch/LM/bioinfo/env/mageck/
TABLE="/data/crispr_sep_23_new/sep23new_sgRNA_count.tsv"
MLE="/data/crispr_sep_23_new/crispr_sep_23_mle_ordered"
P_METHOD="fdr"
DESIGN="/data/crispr_sep_23_new/design_deb_ordered.txt"

mageck mle -k $TABLE -n $MLE --design-matrix $DESIGN --permutation-round 2 --adjust-method fdr --threads 24 --norm-method none
source deactivate
