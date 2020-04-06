import os
import sys

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import re



def load_met(inputfolder, sample, mettype):
    # matches filenames like "footprinting_24_met_cpg.tsv" or "invitro12_2_0_met_dam.tsv"
    filename_regex = re.compile('^(.*)_([0-9]*)_met_([^_]*)\.tsv$')

    ret = None
    for f in os.listdir(inputfolder):
        mitch = filename_regex.match(f)
        if not mitch is None:
            s, batch, m = mitch.groups()
            # Only load select sample and mettypes
            if not s == sample or not m == mettype:
                continue

            print(f)
            met_part = pd.read_csv(os.path.join(inputfolder, f), sep='\t')
            met_part.drop('read_name',axis=1)
            if ret is None:
                ret = met_part.copy()
            else:
                ret = ret.append(met_part)
            break

    return ret

def main(argv):
    inputfolder = argv[1]
    outputfile = argv[2]
    sample = argv[3]
    mettype = argv[4]
    met = load_met(inputfolder, sample, mettype)

    counts,bins  = np.histogram(np.clip(met.log_lik_ratio,-20,20), bins=200)
    metrate = (met.log_lik_ratio > 2.5).sum() / len(met.log_lik_ratio)# (np.abs(met.log_lik_ratio) > 2.5).sum()
    ambiguous = (np.abs(met.log_lik_ratio)<=2.5).sum() / len(met.log_lik_ratio) 
    for i in range(len(counts)):
        loc = (bins[i]+bins[i+1])/2
        c = np.clip(np.abs(loc)/2.5,0,1)
        plt.bar(loc,counts[i],width=bins[i+1]-bins[i], color=[1-c,0.5,c])

    plt.title('LLR %s %s (%d%% methylated, %d%% ambiguous) ' % (sample, mettype, metrate*100, ambiguous*100))
    plt.savefig(outputfile)
   

if __name__ == '__main__':
    main(sys.argv)

