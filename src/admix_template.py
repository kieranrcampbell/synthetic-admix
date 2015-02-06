"""
admix_template.py

Use as a template to create synthetic admix for 2 cell types. Params
<MIXPROP1> : mixing proportion of cell type 1
<MIXPROP2> : 1 - MIXPROP1
<DEPTH> : sequencing depth (e.g. 25000000)
"""
import pandas as pd

dir1 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type1"
dir2 = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/type2"
output_file = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/admix_output/out<STR_MIXPROP>.fastq.gz"

meta_file = "/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/data/meta.csv"
df_meta = pd.read_csv(meta_file,index_col=False)

type_srr_index =  [(df_meta['cell_type'] == i) & (df_meta['usage_type'] == 'bulk') for i in (1,2)]

type_srr = [df_meta['srr_ids'][index].tolist() for index in type_srr_index]



import imp
import numpy as np


createadmix = imp.load_source('createadmix','/net/isi-scratch/kieran/admix/synthetic/synthetic-admix/src/createadmix.py')


sa = createadmix.SyntheticAdmixCreator([dir1, dir2], np.array([<MIXPROP1>, <MIXPROP2>]), 
	<DEPTH>, output_file, type_srr)
sa.cellIO()



