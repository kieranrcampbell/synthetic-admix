# synthetic-admix
Create synthetic bulk RNA-seq data from a mix of single-cell data of known types

## Workflow

The basic workflow is as follows (all scripts can be found in `/src/` directory):

#### 1. Select 40 cells from the `monocle` dataset

The script `monocle_select.R` loads in the monocle dataset. It selects 20 cells from the begining 
of pseudotime (type_1 proliferating cells) and 20 cells from the end of pseudotime (type_2 or differentiating
myoblasts.) It then finds the SRX then SRR IDs of the files and writes two bash scripts, download_type1
and download_type2.

#### 2. Download the corresponding SRA files using `prefetch` from sra-toolkit

Running each bash script will cause prefetch to download the SRA files corresponding to each cell. This uses
(if available) ascp, which is much faster than http. Prefetch downloads to a particular location that can't
be set at runtime (currently `/net/isi-scratch/kieran/ncbi`), so after running these must be moved to
`/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type*`.

#### 3. Convert the SRA files to .fastq.gz using `fastq-dump` from sra-toolkit

The script `fastqmakefilecreate.py` will populate `data/sra_fastqdump_makefiles` with bash scripts to convert
the SRA files to .fastq.gz using `fastq-dump`. It will also create a bash script called makefile to submit each
individual bash script to the cluster, so after running this python script call ./makefile to start the process.

#### 4. Create synthetic .fastq.gz files using `createadmix.py`



#### 5. Compute transcript quantification using `tophat` and `cufflinks`
