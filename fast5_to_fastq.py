import h5py
import sys
import os

def find_fast5(indir, fl=[]):
    ll = os.listdir(indir)
    for l in ll:
        l = os.path.join(indir,l)
        if os.path.isdir(l):
            fl = fl + find_fast5(l, fl)
        elif l.endswith('.fast5'):
            fl.append(l)
    return fl

inputdir=sys.argv[1]
files = find_fast5(inputdir)

for fn in files:
    try:
        with h5py.File(fn,'r') as f:
            for read in f.keys():
                fq = f[read]['Analyses/Basecall_1D_001/BaseCalled_template/Fastq'][()].decode('utf8')
                fq = fq.split('\n')
                if read.startswith('read_'):
                    read = read[5:]
                fq[0] = '@%s'%read

                print('\n'.join(fq), end='')
    except:
        pass
