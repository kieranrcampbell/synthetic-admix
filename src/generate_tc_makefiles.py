"""
Generate bash script makefiles for tophat - cufflinks quantification
"""

import os

TC_bash_template = "tophat -p 32 -o <ALIGNMENT_DIR> -G /net/isi-mirror/ucsc/hg19/hg19_genes.gtf /net/isi-mirror/ucsc/hg19/index/bowtie/hg19_ERCC92 <FASTQ1> <FASTQ2> \n"
TC_bash_template += "cufflinks -o <ABUNDANCE_DIR> -G /net/isi-mirror/ucsc/hg19/hg19_genes.gtf <ALIGNMENT_BAM> \n" 

"""
Need to replace
- <ALIGNMENT_DIR> : which folder should the BAM file be written to?
					(probably singletest/alignment/gamma_x)
- <FASTQX> : X = (1,2) where are the fastq files for this gamma?
					(probably fastq_dirs/mixed_x.fastq for x in (1,2))
- <ALIGNMENT_BAM> : ALIGNMENT_DIR + "/accepted_hits.bam"
- <ABUNDANCE_DIR> : the output directory for this sample

"""

"""
We have:
	- 9 paired end reads in data/admix_output that need quantified
	- 20 paired end reads in data/type1
	- 20 paired end reads in data/type2

"""

admix_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/admix_output"
type1_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type1"
type2_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type2"

all_dirs = [admix_dir, type1_dir, type2_dir]


def get_unique_names(d):
	files = os.listdir(d)
	files = [f for f in files if ".fastq.gz" in f]
	files = [f.replace(".fastq.gz","").split("_")[0] for f in files]
	return list(set(files))

all_files = [get_unique_names(d) for d in all_dirs]


