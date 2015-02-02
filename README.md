# synthetic-admix
Create synthetic bulk RNA-seq data from a mix of single-cell data of known types

## workflow

The basic workflow is as follows:

### 1. Select 40 cells from the `monocle` dataset

### 2. Download the corresponding SRA files using `prefetch` from sra-toolkit

### 3. Convert the SRA files to .fastq.gz using `fastq-dump` from sra-toolkit

### 4. Create synthetic .fastq.gz files using `createadmix.py`

### 5. Compute transcript quantification using `tophat` and `cufflinks`
