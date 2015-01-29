"""
Create makefiles for the cluster for converting sra to fastq

kieran.campbell@dpag.ox.ac.uk
"""

import os

"""
Paths to data - edit to correct
"""
type1_dir = '/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type1'
type2_dir =  '/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type2'
makefile_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/makefiles"


fastq_dump = "/net/isi-software/server/ncbi_sra_toolkit/ncbi_sra_toolkit_2.2.0/bin/fastq-dump" 
opts = "--split-3 --gzip --outdir"

type1_sra_files = [f for f in os.listdir(os.path.join(type1_dir,'sra')) if f.endswith("sra")]
type2_sra_files= [f for f in os.listdir(os.path.join(type2_dir,'sra')) if f.endswith("sra")]

type1_sra = [x.rstrip(".sra") for x in type1_sra_files]
type2_sra = [x.rstrip(".sra") for x in type2_sra_files]

makefile_template = fastq_dump + " " + opts
makefile_template_1 = makefile_template + " " + type1_dir
makefile_template_2 = makefile_template + " " + type2_dir

makefiles = [makefile_template_1 + " " + d for d in [os.path.join(type1_dir,'sra',sra_file) for sra_file in type1_sra_files]]
makefiles = makefiles + [makefile_template_2 + " " + d for d in [os.path.join(type2_dir,'sra',sra_file) for sra_file in type2_sra_files]]

sra = type1_sra + type2_sra
makefile_names = [s + "_makefile" for s in sra]



for i in range(len(makefile_names)):
	pth = os.path.join(makefile_dir, makefile_names[i])
	f = open(pth, 'w')
	f.write(makefiles[i] + "\n")
	f.close()


cluster_template = "nice -19 qsub -q medium_jobs.q -cwd -N "
cluster_cmd = [cluster_template + '"' + s[0] + '" ' + s[1] for s in zip(sra, makefile_names)]

f = open(os.path.join(makefile_dir,"makefile"),'w')
f.write("\n".join(cluster_cmd) + "\n")
f.close()



