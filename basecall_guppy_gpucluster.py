import os
import sys
import shutil
import subprocess

guppy='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/guppy/ont-guppy/bin/guppy_basecaller'
BATCH_SIZE=15
if len(sys.argv) != 3:
    print("Usage: %s inputdir outputdir" % sys.argv[0])
    sys.exit(1)

LSB_JOBID=os.getenv('LSB_JOBID')

inputdir = sys.argv[1]
outputdir = sys.argv[2]

fast5files = [f for f in os.listdir(inputdir) if f.endswith('fast5')]

# Filter out files we already basecalled, so we don't do it twice
if os.path.exists(outputdir):
    done_files = os.listdir(outputdir)
    print('Found %d fast5files. %d files are already basecalled' % (len(fast5files), len(done_files)))
    fast5files = [f for f in fast5files if not f in done_files]
    print('Need to basecall %d files' % len(fast5files))

batches = [fast5files[i:(i+BATCH_SIZE)] for i in range(0,len(fast5files), BATCH_SIZE)]

ssd_in=os.path.join('/ssd/r933r/', LSB_JOBID, 'input')
ssd_out=os.path.join('/ssd/r933r/', LSB_JOBID,'output')

for batch_i in range(len(batches)):
    print("Copying data to SSD (batch %d from %d)" % (batch_i, len(batches)))
    os.mkdir(ssd_in)
    os.mkdir(ssd_out)

    batch = batches[batch_i]

    for f in batch:
        shutil.copyfile(os.path.join(inputdir,f), os.path.join(ssd_in,f))

    guppy_command = '{guppy} --input_path {indir} --save_path {outdir} --device auto ' \
            '--config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg --gpu_runners_per_device 1 ' \
            '--fast5_out --post_out'.format(guppy=guppy, indir=ssd_in, outdir=ssd_out).split(' ')

    p = subprocess.Popen(guppy_command, stdout=subprocess.PIPE)

    for line in p.stdout:
        print('Guppy: ' + line.decode('UTF-8'))

    p.wait()

    # The directory where guppy stores the basecalled fast5 files, since we
    # are only really interested in those and not so much the fastq files
    basecalled_fast5_dir = os.path.join(ssd_out, 'workspace')
    try:
        os.mkdir(outputdir)
    except:
        pass
    print("Copying result to %s" % outputdir)
    for f in os.listdir(basecalled_fast5_dir):
        shutil.copyfile(os.path.join(basecalled_fast5_dir, f), os.path.join(outputdir, f))

    print("Cleaning up %s" % ssd_out)
    shutil.rmtree(ssd_out)
    shutil.rmtree(ssd_in)

with open(os.path.join(outputdir,'done'), 'w') as fp: 
    pass
