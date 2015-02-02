"""

Generate python files + corresponding makefiles for admix pipeline

"""
import numpy as np
import os

makefile_dir = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/synth_fastq_makefiles"


template_file = "admix_template.py"

template = open(template_file,'r')

lines = template.readlines()
# strip whitespace
lines = [l for l in lines if l != '\n']

depth = 25000000
mixprops = np.arange(0.1, 1, 0.1) 
makefile_names = []

for mixprop in mixprops:
	str_mix = str(mixprop).replace(".","")
	fname = os.path.join(makefile_dir,"admix_" + str_mix+ ".py")
	f = open(fname,'w')
	
	for line in lines:
		line = line.replace("<MIXPROP1>",str(mixprop)).replace("<MIXPROP2>",str(1-mixprop))
		line = line.replace("<DEPTH>",str(depth)).replace("<STR_MIXPROP>",str_mix)
		f.write(line)
	f.close()

	fname_makefile = os.path.join(makefile_dir,"makefile" + str_mix)
	makefile_names.append(fname_makefile)
	f = open(fname_makefile,'w')
	f.write("source /net/isi-software/apps/environments/cgat.bash\n")
	f.write("python " + fname + "\n")
	f.close()

f = open(os.path.join(makefile_dir,"makefile"),'w')
for makefile in makefile_names:
	f.write('nice -19 qsub -q medium_jobs.q -cwd -N "' + makefile.split('/')[-1] + '" ' + makefile + '\n')
f.close()

