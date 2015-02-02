"""
testcreateadmix.py
"""

dir1 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type1"
dir2 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type2"
output_file = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/admix_output/out<MIXPROP1>.fastq.gz"

import imp
import numpy as np


createadmix = imp.load_source('createadmix','/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/src/createadmix.py')


sa = createadmix.SyntheticAdmixCreator([dir1, dir2], np.array([0.5, 0.5]), 1000, output_file)
sa.cellIO()



