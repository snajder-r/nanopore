import os
import sys
import pandas as pd
import numpy as np
import re

metcall_dir = sys.argv[1]
sample = sys.argv[2]
chrom = sys.argv[3]
pickle_file = sys.argv[4]

sample_met = None

sample_batch_re = re.compile('(..)_(.*)_.*')


for f in os.listdir(metcall_dir):
    mitch = sample_batch_re.match(f)
    s = mitch.group(1)
    batch = mitch.group(2)
    if s != sample:
        continue
    
    print(f)

    lmet = pd.read_csv(os.path.join(metcall_dir,f), sep='\t', header=0, 
                           dtype={'chromosome':'category', 'strand':'category', 'num_calling_strands':np.uint8, 'num_motifs':np.uint8})

    lmet = lmet.loc[lmet.chromosome==chrom]
    if sample_met is None:
        sample_met = lmet
    else:
        sample_met = pd.concat((sample_met, lmet))

sample_met.to_pickle(pickle_file, compression='gzip')

