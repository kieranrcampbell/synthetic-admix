
# Choose cells for synthetic admix ----------------------------------------
# kieran.campbell@dpag.ox.ac.uk

library(monocle)
library(GEOquery)
library(SRAdb)

set.seed(123)

data(HSMM)

setwd("/net/isi-scratch/kieran/admix/synthetic/cell_selection")

pd <- pData(HSMM)

# Cell states:
# State 1 -> Proliferating cell
# State 2 -> Differentiating myoblast
# State 3 -> Interstitial mesenchmal cell
# 
# Want cells type 1 & 2 so remove type 3 and choose at random
# 20 cells from states 1 & 2

prol_diff <- pd[pd$State != 3,]
pd_sorted <- prol_diff[order(prol_diff$Pseudotime),]

N <- dim(pd_sorted)[1]
state1_cells <- as.character(pd_sorted$cell_id[1:40])
state2_cells <- as.character(pd_sorted$cell_id[(N-39):N])

gse_no <- "GSE52529"
gse <- getGEO(gse_no)

f <- rbind(pData(gse[[1]]), pData(gse[[2]])) 
rownames(f) <- f$description

f1 <- f[state1_cells,]
f2 <- f[state2_cells,]


# use f$description and f$geo_accession

srx_ids <- lapply(list(f1,f2), function(fi) {
  srx_address <- fi$supplementary_file_1
  srx <- sapply(srx_address, function(addr) {
    addr <- as.character(addr)
    split <- strsplit(addr, "/",fixed=T)[[1]][11]
  })
})

sra <- "SRP033135"

if(!file.exists('SRAdb.sqlite')) {
  sqlfile <- getSRAdbFile()
} else {
  sqlfile <- 'SRAmetadb.sqlite'
}

sra_con <- dbConnect('SQLite',sqlfile)
l <- listSRAfile("SRP033135",sra_con)

rownames(l) <- l$experiment
srr_ids <- lapply(srx_ids,function(srx_id) {
  return(l[srx_id,]$run)
})


# Construct metadata file -------------------------------------------------

df_meta <- data.frame(cell_name = c(state1_cells,state2_cells))
df_meta$cell_type <- c(rep(1,40),rep(2,40))
df_meta$cell_type_name <- c(rep("proliferating_cell",40),rep("differentiating_myoblast",40))
df_meta$srx_ids <- unlist(srx_ids)

df_meta$srr_ids <- unlist(srr_ids)

## want to choose half the samples from each type to be used for the synthetic bulk,
## and half to be for single-cell samples. df_meta$usage_type is 'bulk' for those to
## be used as bulk, and 'single' for those to be used as single, decided randomly
## (but with a fixed set seed)

df_meta$usage_type <- "bulk"
single_index <- c(sample(1:40,20),sample(41:80,20))
df_meta$usage_type[single_index] <- "single"



write.csv(df_meta, file="/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/meta.csv",row.names=F)
stop("done")


cmd <- '$PREFETCH -t ascp -a "/home/kieranc/.aspera/connect/bin/ascp|/home/kieranc/.aspera/connect/etc/asperaweb_id_dsa.openssh"'

type1_cmd <- paste(cmd, srr_ids[[1]])
type2_cmd <- paste(cmd, srr_ids[[2]])

# need to go into vdbconfig in sra tools to change download location!
# can be done by <root-to-sra-tools>/bin/vdb-config -i
# here <root-to-sra-tools> = /net/isi-scratch/kieran/tools/sra_toolkit

# set up alias
alias  <- 'PREFETCH=/net/isi-scratch/kieran/tools/sra_toolkit/bin/prefetch'

fileConn1 <- file("download_type1")
writeLines(c(alias,type1_cmd), fileConn1)
close(fileConn1)

fileConn2 <- file("download_type2")
writeLines(c(alias,type2_cmd), fileConn2)
close(fileConn2)

# Then run download_type1, move all to particular directory and repeat for download_type2


## now create makefiles for converting fastq-dumps



