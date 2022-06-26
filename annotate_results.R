# Install & load packages
rm(list = ls())
library(GenomicFeatures)
library(VariantAnnotation)
library(MEDIPS)
library(dplyr)

# Load results from medipAnalysis.R
load(paste(resultsDirectory, comparisonNames, "/methResults.RData", sep=""))

# Make TxDb object

annot <- "GCF_003709585.1_Aci_jub_2_genomic.gff"
TxDb <- makeTxDbFromGFF(file = annot, 
                        dataSource = ann_source, 
                        organism = org_name,
                        taxonomyId = taxID)


# Select DMRs with p-value 0.1
adj_p_val_0.05 <- MEDIPS.selectSig(methResults, p.value = 0.05, adj = T)

# Splitting tables by their methylation status 
# gain = hypermethylation 
gain <- adj_p_val_0.05[which(adj_p_val_0.05$edgeR.logFC > 0),]
# loss = hypomethylation
loss <- adj_p_val_0.05[which(adj_p_val_0.05$edgeR.logFC < 0),]

## Formatting the data for the annotations
gain_m <- MEDIPS.mergeFrames(gain)
loss_m <- MEDIPS.mergeFrames(loss)

gain_mg <- makeGRangesFromDataFrame(gain_m, keep.extra.columns = T)
loss_mg <- makeGRangesFromDataFrame(loss_m, keep.extra.columns = T)

TxDb_gain <- TxDb_loss  <- TxDb 

## Annotation
# Hypermethylated tables p-val 0.05
seqlevels(TxDb_gain) <- intersect(seqlevels0(TxDb), seqlevels(gain_mg))

gain_ann <- as_tibble(locateVariants(gain_mg, TxDb_gain, AllVariants())) %>% 
  select_if(~ !is.list(.)) %>% 
  mutate(QUERYID = paste0(QUERYID, "_hyper"))


# Hypomethylated tables p-val 0.05
seqlevels(TxDb_loss) <- intersect(seqlevels0(TxDb), seqlevels(loss_mg))

loss_ann <- as_tibble(locateVariants(loss_mg, TxDb_loss, AllVariants())) %>% 
  select_if(~ !is.list(.)) %>%
  mutate(QUERYID = paste0(QUERYID, "_hypo"))


# Create Table for hypo and hypermethylated regions
annotated_result_p.005 <- bind_rows(as_tibble(gain_ann), as_tibble(loss_ann))

write.csv2(annotated_result_p.005, file = "annot_005.csv")


#####same with p = 05

# rm(list = ls())
# # Load results from medipAnalysis.R
# load("/data/fg2/vullioud/genetic_hyena/results/EdgeR_fit_final_300.rda")
# 
# # Make TxDb object
# #annot <- "/data/fg2/weyrich/Hyena/2021_MBD_Seq/MEDIPS/data/refGenome/Genome_Annotation_spottedHyena_UniP/Hyena_Liftoff_mapAnnotation_ref_Hhy_NCBI_ASM300989v1/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
# annot <- "/data/shared/AlexWeyrich/refGenome_medips/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
# TxDb <- makeTxDbFromGFF(file = annot, 
#                         dataSource = "Hofreiter_UniP", 
#                         organism = "Crocuta crocuta")
# 
# 
# # Select DMRs with p-value 0.1
# adj_p_val_0.05 <- MEDIPS.selectSig(EdgeR_fit_final_300$result, p.value = 0.05, adj = T)
# 
# # Splitting tables by their methylation status 
# # gain = hypermethylation 
# gain <- adj_p_val_0.05[which(adj_p_val_0.05$edgeR.logFC > 0),]
# # loss = hypomethylation
# loss <- adj_p_val_0.05[which(adj_p_val_0.05$edgeR.logFC < 0),]
# 
# ## Formatting the data for the annotations
# gain_m <- MEDIPS.mergeFrames(gain)
# loss_m <- MEDIPS.mergeFrames(loss)
# 
# gain_mg <- makeGRangesFromDataFrame(gain_m, keep.extra.columns = T)
# loss_mg <- makeGRangesFromDataFrame(loss_m, keep.extra.columns = T)
# 
# TxDb_gain <- TxDb_loss  <- TxDb #back up to intersect the seqlevels in next step 
# ## Annotation
# # Hypermethylated tables p-val 0.05
# seqlevels(TxDb_gain) <- intersect(seqlevels0(TxDb), seqlevels(gain_mg))
# 
# gain_ann <- as_tibble(locateVariants(gain_mg, TxDb_gain, AllVariants())) %>% 
#   select_if(~ !is.list(.)) %>% 
#   mutate(QUERYID = paste0(QUERYID, "_hyper"))
# 
# 
# # Hypomethylated tables p-val 0.05
# seqlevels(TxDb_loss) <- intersect(seqlevels0(TxDb), seqlevels(loss_mg))
# 
# loss_ann <- as_tibble(locateVariants(loss_mg, TxDb_loss, AllVariants())) %>% 
#   select_if(~ !is.list(.)) %>%
#   mutate(QUERYID = paste0(QUERYID, "_hypo"))
# 
# 
# # Create Table for hypo and hypermethylated regions11
# annotated_result_p.05 <- bind_rows(as_tibble(gain_ann), as_tibble(loss_ann))
# 
# save(annotated_result_p.05, file = "/data/fg2/vullioud/genetic_hyena/results/annot_edgeR_300_p05.rda")

