#import libraries
library(CRISPRcleanR)
library(tidyverse)

#Load libraries
data("KY_Library_v1.0")
#'%!in%' <- function(x,y)!('%in%'(x,y))
fn = "/data/crispr_sep_23_new/crispr_sep_23.count.txt"

#Read CRISPR counts
df <- read.delim(fn, sep = "\t", header = T)

#Normalization of sgRNA counts and fold change computation
##Adjust for library size and read count distributions, scaling by the tot number of reads peer sample 
normANDfcs <- ccr.NormfoldChanges(fn,
                                  min_reads = 15,           
                                  EXPname = 'crispr_sep_23',
                                  libraryAnnotation = KY_Library_v1.0,
                                  saveToFig = T,
                                  outdir = "sep23/",
                                  ncontrols = 3,             #number of replicates 
                                  method='ScalingByTotalReads')

##sorting sgRNA according chromosomes 
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,KY_Library_v1.0)

##Correction for gene independent responses to CRISPR-Cas9 targeting
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label='crispr_sep_23',saveTO = "sep23/")

#Visualization and assessment of results
data("BAGEL_essential")
data("BAGEL_nonEssential")

##Store logFCs
FCs <- correctedFCs$corrected_logFCs$avgFC
names(FCs) <- rownames(correctedFCs$corrected_logFCs)

##Convert reference sets of CFE(core-fitness essential) and non-essential genes into sets of sgRNAs
BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_nonEssential)

##ROC curve quantifying the performances in classifying the considered reference sets based on their logFCs
sgRNA_level_ROC <- ccr.ROC_Curve(FCs, BAGEL_essential_sgRNAs,BAGEL_nonEssential_sgRNAs)
sgRNA_level_ROC$curve

##Gene level analysis, for each gene average logFC across its targeting sgRNAs
geneFCs <- ccr.geneMeanFCs(FCs, KY_Library_v1.0)
gene_level_ROC<-ccr.ROC_Curve(geneFCs, BAGEL_essential, BAGEL_nonEssential,FDRth = 0.05)

##Precision/Recall evaluation, quantification of performances in classifying the considered reference sets based in their logFC
sgRNA_level_PrRc <- ccr.PrRc_Curve(FCs, BAGEL_essential_sgRNAs, BAGEL_nonEssential_sgRNAs)

##Recall values at fixed False Discovery Rate 
gene_level_Pr_Rc <- ccr.PrRc_Curve(geneFCs, BAGEL_essential,BAGEL_nonEssential,FDRth = 0.05)


#Depletion profile visualization with genes signature superimposed and recall computation
data("EssGenes.ribosomalProteins")
data("EssGenes.DNA_REPLICATION_cons")
data("EssGenes.KEGG_rna_polymerase")
data("EssGenes.PROTEASOME_cons")
data("EssGenes.SPLICEOSOME_cons")

SIGNATURES <- list(Ribosomal_Proteins = EssGenes.ribosomalProteins,
                   DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                   RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                   Proteasome = EssGenes.PROTEASOME_cons,
                   Spliceosome = EssGenes.SPLICEOSOME_cons,
                   CFE = BAGEL_essential,
                   non_essential = BAGEL_nonEssential)

##Visualization of the gene essentiality profile
recall_scores <- ccr.VisDepAndSig(FCsprofile = geneFCs,
                                  SIGNATURES = SIGNATURES,
                                  TITLE = 'sep23',
                                  pIs =6,
                                  nIs = 7,
                                  th = 0.05)

#Recall variations following CRISPRcleanR correction for reference, copy number amplified, and non expressed genes
ccr.RecallCurves(correctedFCs$corrected_logFCs, libraryAnnotation = KY_Library_v1.0, cellLine = "A549")


#Save normalized counts
colnames(normANDfcs$norm_counts)[2] <- "gene"

correctedCounts <- ccr.correctCounts('sep23',
                                     normANDfcs$norm_counts,
                                     correctedFCs,
                                     KY_Library_v1.0,
                                     minTargetedGenes = 3,
                                     OutDir = 'sep23/',
                                     ncontrols = 3)

ccr.PlainTsvFile(correctedCounts, "sep23", "sep23/")


