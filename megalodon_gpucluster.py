import os
import sys
import shutil
import subprocess

guppy='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/guppy/ont-guppy/bin/guppy_basecall_server'
megalodon='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/bin/megalodon'
calibration='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/lib/python3.7/site-packages/megalodon/model_data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg/megalodon_mod_calibration.npz'

BATCH_SIZE=1
if len(sys.argv) != 3:
    print("Usage: %s inputdir outputdir" % sys.argv[0])
    sys.exit(1)

LSB_JOBID=os.getenv('LSB_JOBID')

inputdir = sys.argv[1]
outputdir = sys.argv[2]

fast5files = [f for f in os.listdir(inputdir) if f.endswith('fast5')]

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

    megalodon_command = '{megalodon} --overwrite --devices 0 --outputs=mod_basecalls --guppy-server-path={guppy} --output-directory={outdir} ' \
            '--mod-calibration-filename={calibration} {indir}'.format(
                    guppy=guppy, megalodon=megalodon, indir=ssd_in, 
                    outdir=ssd_out, calibration=calibration).split(' ')

    print(megalodon_command)

    p = subprocess.Popen(megalodon_command, stdout=subprocess.PIPE)

    for line in p.stdout:
        print('Guppy: ' + line.decode('UTF-8'))

    p.wait()

    print("Copying result to %s" % outputdir)
    try:
        os.mkdir(outputdir)
    except:
        pass

    shutil.copytree(ssd_out, os.path.join(outputdir,str(batch_i)))
    #shutil.copy(os.path.join(ssd_out,'guppy_log.out'), os.path.join(outputdir,'guppy_%d.log'%batch_i))
    #shutil.copy(os.path.join(ssd_out,'guppy_log.err'), os.path.join(outputdir,'guppy_%d.err'%batch_i))
    #shutil.copy(os.path.join(ssd_out,'basecalls.modified_base_scores.hdf5'), os.path.join(outputdir,'modified_base_scores_%d.hdf5'%batch_i))

    print("Cleaning up %s" % ssd_out)
    shutil.rmtree(ssd_out)
    shutil.rmtree(ssd_in)
