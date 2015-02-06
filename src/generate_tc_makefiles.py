"""
Generate bash script makefiles for tophat - cufflinks quantification
"""

import os
import pandas as pd

TC_bash_template = "export PATH=/net/isi-scratch/kieran/tools/bowtie2-2.2.4:$PATH\n"
TC_bash_template += "TOPHAT=/net/isi-scratch/kieran/tools/tophat-2.0.13.Linux_x86_64/tophat\n"
TC_bash_template += "CUFFQUANT=/net/isi-scratch/kieran/tools/cufflinks-2.2.1.Linux_x86_64/cuffquant\n"
TC_bash_template += "$TOPHAT -p 4 -o <ALIGNMENT_DIR> -G /net/isi-mirror/ucsc/hg19/hg19_genes.gtf /net/isi-mirror/ucsc/hg19/index/bowtie/hg19_ERCC92 <FASTQ1> <FASTQ2> \n"
TC_bash_template += "$CUFFQUANT -o <ABUNDANCE_DIR> -p 4 /net/isi-mirror/ucsc/hg19/hg19_genes.gtf <ALIGNMENT_BAM> \n" 

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

meta_file = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/meta.csv"
df_meta = pd.read_csv(meta_file,index_col=False)

""" Only want to compute expression estimates for single cells not included in bulk estimates """
single_srr =  df_meta['srr_ids'][df_meta['usage_type'] == 'single'].tolist()


data_path = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/"
cell_types = ["admix_output","type1","type2"]

all_dirs = [os.path.join(data_path, c) for c in cell_types]


def get_unique_names(d):
	files = os.listdir(d)
	files = [f for f in files if ".fastq.gz" in f]
	files = [f.replace(".fastq.gz","").split("_")[0] for f in files]
	return list(set(files))

all_files = [get_unique_names(d) for d in all_dirs]

""" strip out any bulk examples """
all_files[1] = list(set(all_files[1]).intersection(set(single_srr)))
all_files[2] = list(set(all_files[2]).intersection(set(single_srr)))


align_dir = os.path.join(data_path, "transcript_quant","alignments")
abund_dir = os.path.join(data_path, "transcript_quant","abundance")
align_dirs = []
abund_dirs = []

assert len(all_files) == len(all_dirs), "Length of directories must match"
for i in range(len(all_files)):
	align_dirs.append([])
	abund_dirs.append([])
	for f in all_files[i]:
		full_path_align = os.path.join(align_dir, cell_types[i],f)
		if not os.path.exists(full_path_align):
			os.makedirs(full_path_align)
		align_dirs[i].append(full_path_align)

		full_path_abund = os.path.join(abund_dir, cell_types[i],f)
		if not os.path.exists(full_path_abund):
			os.makedirs(full_path_abund)
		abund_dirs[i].append(full_path_abund)


fastq_paths = [[os.path.join(all_dirs[i],f) for f in all_files[i]] for i in range(len(all_dirs))]

all_files_flat = [f for path in all_files for f in path]
fastq_paths_flat = [f for path in fastq_paths for f in path]
align_dirs_flat = [f for path in align_dirs for f in path]
abund_dirs_flat = [f for path in abund_dirs for f in path]

all_paths = zip(fastq_paths_flat, align_dirs_flat, abund_dirs_flat)

makefiles = {}
for i in range(len(all_paths)):
	paths = all_paths[i] # fastq - align - abund
	mkfile = TC_bash_template.replace("<ALIGNMENT_DIR>",paths[1])
	mkfile = mkfile.replace("<ABUNDANCE_DIR>",paths[2])
	mkfile = mkfile.replace("<ALIGNMENT_BAM>",os.path.join(paths[1],"accepted_hits.bam"))
	fastqs = [paths[0] + "_" + str(j) + ".fastq.gz" for j in (1,2)]
	mkfile = mkfile.replace("<FASTQ1>",fastqs[0])
	mkfile = mkfile.replace("<FASTQ2>",fastqs[1])
	makefiles[all_files_flat[i] + "_makefile"] = mkfile

makefile_dir = os.path.join(data_path,"transcript_quant","makefiles")
for makefile, content in makefiles.iteritems():
	f = open(os.path.join(makefile_dir, makefile),'w')
	f.writelines(content)
	f.close()

f = open(os.path.join(makefile_dir, "makefile"),'w')

for makefile in makefiles:
	process_name = makefile.replace("_makefile","")[-4:]
	process_name = "top_cuff_" + process_name
	f.write('nice -19 qsub -pe dedicated 4 -q medium_jobs.q -cwd -N "' + process_name + '" ' + makefile + "\n")
f.close()

bulk_abund = abund_dirs[0]
single_abund = abund_dirs[1] + abund_dirs[2]
bulk_names = all_files[0]
single_names = all_files[1] + all_files[2]


bulk = [os.path.join(s,"abundances.cxb") for s in bulk_abund]
single = [os.path.join(s,"abundances.cxb") for s in single_abund]

sample_sheet_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/transcript_quant/norm_abund"

def write_sample_sheet(file_name, cxb_paths, cell_names):
	f = open(os.path.join(sample_sheet_dir, file_name),'w')
	f.write("sample_name\tgroup\n")
	for item in zip(cxb_paths,cell_names):
		f.write(item[0] + "\t" + item[1] + "\n")
	f.close()

fnames = ["bulk_sample_sheet.txt","single_sample_sheet.txt"]

write_sample_sheet(fnames[0], bulk, bulk_names)
write_sample_sheet(fnames[1], single, single_names)

outdirs = ["bulk","single"]

for d in outdirs:
	if not os.path.join(sample_sheet_dir, d):
		os.makedirs(os.path.join(sample_sheet_dir,d))

		

cmd = "CUFFNORM=/net/isi-scratch/kieran/tools/cufflinks-2.2.1.Linux_x86_64/cuffnorm\n"
cmd += "$CUFFNORM --use-sample-sheet -o <outdir> /net/isi-mirror/ucsc/hg19/hg19_genes.gtf <sample_sheet>\n"

makefile_names = []

for i in (0,1):
	c = cmd.replace("<outdir>",os.path.join(sample_sheet_dir,outdirs[i]))
	c = c.replace("<sample_sheet>",os.path.join(sample_sheet_dir,fnames[i]))
	makefile_name = os.path.join(sample_sheet_dir,outdirs[i] + "_makefile")
	makefile_names.append(makefile_name)
	f = open(makefile_name,'w')
	f.write(c)
	f.close()

f = open(os.path.join(sample_sheet_dir,"makefile"),'w')
for n in makefile_names:
	f.write('nice -19 qsub -q medium_jobs.q -cwd -N "cuffnorm" '+ n + "\n")	
f.close()


