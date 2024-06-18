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
run_dir1="/data/FASTQ/220720_A00302_0433_BHC5Y2DMXY/"
run_dir2="/data/220801_A00302_0435_BHHCKNDSX3/"
FASTQ="$run_dir1/Sample_S43704_T0_1/S43704_T0_1_S78_L001_R1_001.fastq.gz $run_dir1/Sample_S43705_T0_2/S43705_T0_2_S24_L001_R1_001.fastq.gz $run_dir1/Sample_S43706_T0_3/S43706_T0_3_S118_L001_R1_001.fastq.gz $run_dir1/Sample_S43707_14_gg_1/S43707_14_gg_1_S100_L001_R1_001.fastq.gz $run_dir1/Sample_S43708_14gg_2/S43708_14gg_2_S36_L002_R1_001.fastq.gz $run_dir1/Sample_S43709_14_gg_3/S43709_14_gg_3_S19_L002_R1_001.fastq.gz $run_dir1/Sample_S43711_IFN_2/S43711_IFN_2_S80_L002_R1_001.fastq.gz $run_dir1/Sample_S43712_IFN_3/S43712_IFN_3_S63_L002_R1_001.fastq.gz $run_dir2/Sample_S43710_IFN_1/S43710_IFN_1_S2_L002_R1_001.fastq.gz $run_dir2/Sample_S43716_21gg_1/S43716_21gg_1_S98_L001_R1_001.fastq.gz $run_dir2/Sample_S43717_21gg_2/S43717_21gg_2_S177_L001_R1_001.fastq.gz $run_dir2/Sample_S43713_COMBO_1/S43713_COMBO_1_S329_L002_R1_001.fastq.gz $run_dir2/Sample_S43714_COMBO_2/S43714_COMBO_2_S73_L001_R1_001.fastq.gz $run_dir2/Sample_S43715_COMBO_3/S43715_COMBO_3_S314_L001_R1_001.fastq.gz"
out_dir="/data/crispr_sep_23_new/"
OUTPUT_PREFIX="$out_dir/crispr_sep_23"
BASES_TO_TRIM="23"
SGRNA_LIST="/data/ky_lib.csv"
SAMPLE_LABELS="S43704_T0_1,S43705_T0_2,S43706_T0_3,S43707_14_gg_1,S43708_14gg_2,S43709_14_gg_3,S43711_IFN_2,S43712_IFN_3,S43710_IFN_1,S43716_21gg_1,S43717_21gg_2,S43713_COMBO_1,S43714_COMBO_2,S43715_COMBO_3"


mageck count -l $SGRNA_LIST --fastq $FASTQ -n $OUTPUT_PREFIX --trim-5 $BASES_TO_TRIM --sample-label $SAMPLE_LABELS 

source deactivate 
