
# Take the read counts reported by cuffnorm, convert them to TPM
# and output them to count tables

single_file  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/transcript_quant/norm_abund/single/genes.fpkm_table"
bulk_file  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/transcript_quant/norm_abund/bulk/genes.fpkm_table"
meta_file  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/meta.csv"

single_data  <- read.delim(single_file, row.names=1)
bulk_data  <- read.delim(bulk_file, row.names=1)

df_meta <- read.csv(meta_file)

single_names <- sapply(colnames(single_data),function(x) strsplit(x,"_")[[1]][1],USE.NAMES=F)
colnames(single_data) <- single_names

# Convert FPKM -> TPM -----------------------------------------------------

# To convert FPKM to TPM we use
# 
# TPM_i = FPKM_i * 10^6 / \sum_j FPKM_j 
# 
# for a given sample. See 
# https://liorpachter.wordpress.com/2014/04/30/estimating-number-of-transcripts-from-rna-seq-measurements-and-why-i-believe-in-paywall/#comments

fpkm_to_tpm <- function(fpkm) round(fpkm * 1000000 / sum(fpkm))

single_data <- apply(single_data, 2, fpkm_to_tpm)
bulk_data <- apply(bulk_data, 2, fpkm_to_tpm)

type1_names <- as.character(df_meta$srr_ids[df_meta$cell_type == 1 & df_meta$usage_type == "single"])
type2_names <- as.character(df_meta$srr_ids[df_meta$cell_type == 2 & df_meta$usage_type == "single"])

type1_data <- single_data[,type1_names]
type2_data <- single_data[,type2_names]

bulk_output  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/counts/bulk.csv"
type1_output  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/counts/type1.csv"
type2_output  <- "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/counts/type2.csv"

write.csv(bulk_data, bulk_output)
write.csv(type1_data, type1_output)
write.csv(type2_data, type2_output)

