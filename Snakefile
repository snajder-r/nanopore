import os
import sys

configfile: "medulloblastoma.yaml"

for k in config.keys():
    globals()[k] = config[k]


samples = []
batches = []

samplebatches = dict()
def update_batches():
    for s in justsamples:
        samplebatches[s] = []
        batchdir=os.path.join(basedir,'raw', s,'batched')
        if os.path.exists(batchdir):
            for b in os.listdir(batchdir):
                if os.path.isdir(os.path.join(batchdir,b)):
                    samples.append(s)
                    batches.append(b)
                    samplebatches[s].append(b)

def find_fast5(indir, fl=[]):
    ll = os.listdir(indir)
    for l in ll:
        l = os.path.join(indir,l)
        if os.path.isdir(l):
            fl = fl + find_fast5(l, fl)
        elif l.endswith('.fast5'):
            fl.append(l)
    return fl


update_batches()
print(samplebatches)

rule split_batches:
    output: directory(os.path.join(basedir,'raw/{sample}/batched/'))
    run: 
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        i = 0
        b = 0
        sample_dir = os.path.join(basedir,'raw/%s/guppy/' % wildcards.sample)
        print(sample_dir)
        raw_batches = os.listdir(sample_dir)
        while i < len(raw_batches):
            batchdir=os.path.join(output[0], '%d'%b)
            os.mkdir(batchdir)
            for _ in range(per_batch):
                if i == len(raw_batches):
                     break
                src=os.path.join(basedir, os.path.join(sample_dir, raw_batches[i]))
                dst=os.path.join(batchdir, '%s' % raw_batches[i])
                os.symlink(src,dst)
                i+=1
            b+=1
        update_batches()

rule all_split_batches:
    input: expand(rules.split_batches.output, sample=justsamples)

rule fast5_to_fastq:
    input: os.path.join(basedir, 'raw/{sample}/batched/{batch}/')
    output: os.path.join(basedir, 'fastq/{sample}_{batch}.fq')
    params:
        jobname='fq_{sample}{batch}',
        runtime='03:00',
        memusage='8000',
        slots='1',
        misc=''
    run:
        files = find_fast5(input[0])
        with file(output[0],'w') as of:
            for fn in files:
                try:
                    with h5py.File(fn,'r') as f:
                        for read in f.keys():
                            fq = f[read]['Analyses/{basecall_id}/BaseCalled_template/Fastq'.format(basecall_id=basecall_id)][()].decode('utf8')
                            fq = fq.split('\n')
                            if read.startswith('read_'):
                                read = read[5:]
                            fq[0] = '@%s'%read

                            of.write('\n'.join(fq), end='\n')
                except:
                    print('WARN: Could not read basecalls in %s'%fn)
                    pass

rule alignment:
    input: rules.fast5_to_fastq.output
    output: os.path.join(basedir, 'mapping/{sample}_{batch}.sorted.filtered.bam')
    params:
        jobname='align_{sample}_{batch}',
        runtime='09:00',
        memusage='2000',
        slots='1',
        misc=''
    shell: '''
           module load samtools
           minimap2 -a -x map-ont {reference} {input} | samtools view -q 1 -b | samtools sort -T {output}.tmp -o {output}
           samtools index {output}
           '''
rule nanopolish_index:
    input:
        f5=os.path.join(basedir, 'raw/{sample}/batched/{batch}/'),
        fq=rules.fast5_to_fastq.output
    output: os.path.join(basedir, 'fastq/{sample}_{batch}.fq.index.readdb')
    params:
        jobname='npidx_{sample}_{batch}',
        runtime='08:00',
        memusage='2000',
        slots='1',
        misc=''
    shell: 'nanopolish index -d {input.f5} {input.fq}'

rule metcall:
    input: 
        fq=rules.fast5_to_fastq.output,
        bam=rules.alignment.output,
        fqidx=rules.nanopolish_index.output
    output: os.path.join(basedir, 'met/{sample}_{batch}_met_{mtype}.tsv')
    params:
        jobname='metcall_{sample}_{batch}',
        runtime='08:00',
        memusage='8000',
        slots='8 -R "span[hosts=1]"',
        misc=''
    shell: 'nanopolish call-methylation -t 8 -g {reference} -b {input.bam} -r {input.fq} -q {wildcards.mtype} > {output}'

def prepare_mergebams_input(wildcards):
    return expand((os.path.join(basedir, 'mapping/%s' % wildcards.sample) + '_{batch}.sorted.filtered.bam'), batch=samplebatches[wildcards.sample])

rule prepare_mergebams:
    input: prepare_mergebams_input
    output: os.path.join(basedir, 'mapping/{sample}.bam.filelist.txt')
    shell: 'echo {input} | sed \'s/ /\\n/g\' > {output}'

rule mergebams:
    input: rules.prepare_mergebams.output
    output: os.path.join(basedir, 'mapping/{sample}.sorted.bam')
    params:
        jobname='mergebam_{sample}',
        runtime='09:59',
        memusage='32000',
        slots='1',
        misc=''
    shell: '''
           module load samtools
           samtools merge -b {input} - | samtools calmd -b - {reference} | samtools sort -T {output}.tmp -o {output}
           '''

rule varcall:
    input:
        bam=rules.mergebams.output
    output: os.path.join(basedir, 'var/{sample}_sv.vcf')
    params:
        jobname='varcall_{sample}',
        runtime='24:00',
        memusage='16000',
        slots='16 -R "span[hosts=1]"',
        misc=''
    shell: '''
           module load gcc/5.4.0
           sniffles -s 4 -t 16 -m {input} -v {output}
           '''

rule megalodon:
    input: os.path.join(basedir, 'raw/{sample}/batched/{batch}/')
    output: os.path.join(basedir, 'megalodon/{sample}/{batch}/basecalls.modified_base_scores.hdf5')
    params: 
        jobname='megalodon_{sample}_{batch}',
        runtime='48:00',
        memusage='16000',
        slots='16 -R "span[hosts=1]"',
        misc=''
    shell: '{megalodon} --guppy-params "--cpu_threads_per_caller 16" --guppy-server-path {guppy} --overwrite --output-directory %s --mod-calibration-filename {megalodon_calibration} --outputs=mod_basecalls {input}' % os.path.join(basedir, 'megalodon/{wildcards.sample}/{wildcards.batch}/')

def merge_met_perchrom_input(wildcards):
    return expand(os.path.join(basedir, 'met/%s_{batch}_met_%s.tsv' % (wildcards.sample, wildcards.mettype)), batch=samplebatches[wildcards.sample])

rule merge_met_perchrom:
    input: merge_met_perchrom_input
    output: os.path.join(basedir, 'met_merged/{sample}_chr{chrom}_met_{mettype}.pkl')
    params:
        jobname='mergemet_{chrom}',
        runtime='8:00',
        memusage='32000',
        slots='1',
        misc=''
    run:
        metcall_dir = os.path.join(basedir, 'met') 
        sample = wildcards['sample']
        chrom = wildcards['chrom']
        mettype = wildcards['mettype']
        pickle_file = output[0]

        sample_met = None

        sample_batch_re = re.compile('(..)_(.*)_met_(.*)\.tsv')

        for f in os.listdir(metcall_dir):
            mitch = sample_batch_re.match(f)
            s = mitch.group(1)
            batch = mitch.group(2)
            mt = mitch.group(3)
            if s != sample or mettype != mt:
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

rule tombo_make_index:
    input: os.path.join(basedir, 'raw/{sample}/single/')
    output: os.path.join(basedir, 'raw/{sample}/.single.RawGenomeCorrected_000.tombo.index')
    params:
        jobname='tombo_make_index',
        runtime='24:00',
        memusage='16000',
        slots='1',
        misc=''
    shell: '{tombo} filter clear_filters --fast5-basedirs {input}'

rule tombo_resquiggle:
    input: 
        index=rules.tombo_make_index.output,
        directory=rules.tombo_make_index.input
    output: os.path.join(basedir, 'raw/{sample}/.single.resquiggled')
    params:
        jobname='tombo_resquiggle',
        runtime='48:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]',
        misc=''
    shell: '''
           {tombo} resquiggle --processes 32 --dna {input.directory} {reference};
           RESULT=$?
           if $RESULT -eq 0; then touch {output}; fi
           exit $RESULT
           '''

rule tombo_metcall:
    input: resquiggled=rules.tombo_resquiggle.output,
           directory=rules.tombo_make_index.input
    output: stats_6ma=os.path.join(basedir, 'tombo/{sample}.6mA.tombo.stats'),
            stats_5mc=os.path.join(basedir, 'tombo/{sample}.5mC.tombo.stats')
    params:
        jobname='tombo_metcall',
        runtime='48:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]',
        misc=''
    shell: '''
           {tombo} detect_modification alternative_model --fast5-basedirs {input.directory} --statistics-file-basename {wildcards.sample} --alternate-bases 6mA 5mC --processes 32
           '''  



def allmetcall_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        # All mtypes
        for j in range(len(args[2])):
            #      sample   -  batch   -  mtype
            yield args[0][i],args[1][i],args[2][j]
        
def report_methylation_input(wildcards):
    return expand(rules.metcall.output, batch=samplebatches[wildcards.sample], sample=wildcards.sample, mtype=wildcards.mtype)

rule report_methylation:
    input: report_methylation_input
    output: os.path.join(basedir, 'report/{sample}_{mtype}.pdf')
    params:
        jobname='report_{sample}_{mtype}',
        runtime='0:10',
        memusage='4000',
        slots='1',
        misc=''
    shell: '{python} report_met.py ' + os.path.join(basedir,'met') + ' {output} {wildcards.sample} {wildcards.mtype}'

rule report_all_methylation:
    input: expand(rules.report_methylation.output, mtype=['cpg'], sample=justsamples)


rule allmetcall:
    input: expand(os.path.join(basedir, 'met/{sample}_{batch}_met_{mtype}.tsv'), \
               allmetcall_combinator, sample=samples, batch=batches, mtype=['cpg'])

def allmegalodon_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        #      sample   -  batch   -  mtype
        yield args[0][i],args[1][i]


rule allmegalodon:
    input: expand(os.path.join(basedir, 'megalodon/{sample}/{batch}/basecalls.modified_base_scores.hdf5'), \
               allmegalodon_combinator, sample=samples, batch=batches)

rule allmetmerged:
    input: expand(os.path.join(basedir, 'met_merged/{sample}_chr{chrom}_met_cpg.pkl'), sample=justsamples, chrom=chroms)


rule allvarcall:
    input: expand(os.path.join(basedir, 'var/{sample}_sv.vcf'), sample=justsamples)

rule all_tombo_metcall:
    input: expand(rules.tombo_metcall.output.stats_5mc, sample=justsamples)

rule test:
    input: os.path.join(basedir, 'megalodon/invitro11_1/0/basecalls.modified_base_scores.hdf5')

localrules: prepare_mergebams, split_batches
