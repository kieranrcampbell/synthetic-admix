"""
testcreateadmix.py
"""

dir1 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type1"
dir2 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type2"

import createadmix
import numpy as np

reload(createadmix)

sa = createadmix.SyntheticAdmixCreator([dir1, dir2], np.array([0.5, 0.5]), 1000, "output.fastq.gz")



