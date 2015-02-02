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

The python module `createadmix.py` can take a set of fastq.gz files corresponding to different cell types
(each cell type must be a different directory, but there is no limit on the number of cell types supported) and will
create synthetic bulk RNA-seq samples based on mixing the cell types in known ratios. To use this, run `generate_admix_makefiles.py` which will create python scripts based on `admix_template.py` to call `createadmix.py` 
on different mixture ratios (currently 0.1, 0.2,...,0.9). It will also write bash scripts to call these python scripts and a final bash script `makefile` that submits all jobs to the cluster. The directory for all generated python and bash scripts is `/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/synth_fastq_makefiles`.


#### 5. Compute transcript quantification using `tophat` and `cufflinks`

TBC

### Scripts associated with each stage

Stage | Important scripts | Target directory
-------| ---------------------| ------------------
Cell selecting | `monocle_select.R` | Currently `src/`
Downloading SRA | `download_typeA`, `download_typeB` | `/net/isi-scratch/kieran/ncbi`, `data/type1`,`data/type2`
Converting SRA to fastq | `fastqmakefilecreate.py` | `data/sra_fastqdump_makefiles`
Creating synthetic data | `createadmix.py`, `admix_template.py`, `generate_admix_makefiles.py` | `data/synth_fastq_makefiles`, `/data/admix_output`

