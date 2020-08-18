import os
from string import Formatter

for k in config.keys():
    globals()[k] = config[k]

# Mandatory globals
basecall_id = config['basecall_id']
basedir = config['basedir']
mettypes = config['mettypes']
per_batch = config['per_batch']



basecall_group = 'Basecall_1D_%s' % basecall_id
if 'unique_samples' not in config.keys():
    unique_samples = os.listdir(os.path.join(basedir, 'raw'))
if 'fastq_ending' not in config.keys():
    fastq_ending = 'fq'
if 'chroms' not in config.keys():
    chroms = []
if 'is_transcriptome' not in config.keys():
    is_transcriptome = False

'''
Detect batches for each sample. Fills the global variables:
  samples: flat array of length num_samples x num_batches
  batches: flat array of length num_samples x num_batches
  samplebatches: dict containing the list of batches per sample

samples and batches are filled such that zip(samples,batches) would result
in the full list of sample and batch tuples.
'''

sb = glob_wildcards(os.path.join(basedir, 'raw', '{sample}', '{batch}.%s'%fastq_ending))
sbf = glob_wildcards(os.path.join(basedir, 'raw', '{sample}', 'batched', '{batch}', '{filename}.fast5'))

def samplebatches(sample):
    if len(sb.sample) == 0:
        return list(set([sbf.batch[i] for i in range(len(sbf.batch)) if sbf.sample[i] == sample]))
    else:
        return list(set([sb.batch[i] for i in range(len(sb.batch)) if sb.sample[i] == sample]))


'''
##############################################################################
# Combinators 
##############################################################################

Below we define some combinators of wildcards. These can be used in "all_"
type rules. The expand function in Snakemake assumes that the list of 
wildcard combinations is always the full cross-product. In our case, this
is not the case, as samples have different number of batches. 

Therefore we define a few custom combinators below
'''

''' 
This combinator assumes that all wildcard variables are simply zipped 
together. This is designed for the sample and batch variables.
'''


def zip_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        yield tuple(args[j][i] for j in range(len(args)))


'''
This combinator zips the first two arguments, and combines those with the 
third argument. We use this to create all sample-batch-mtype combinations.
'''


def zip2_comb3_combinator(*args, **kwargs):
    # First two wildcards
    for i in range(len(args[0])):
        # Combine with third wildcard
        for j in range(len(args[2])):
            #     wc1   -     wc2   -     wc3
            yield args[0][i], args[1][i], args[2][j]


def globspand(path, **kwargs):
    """
    This is a handy utility function that combines the functionality of expand
    and glob_wildcards. As an input, it takes a path containing wildcards.
    It returns a curried function that takes a wildcards object as an input.
    Unlike expand, however, it can deal with incomplete wildcards. It will then
    complete all missing wildcards using the glob_wildcards function.

    Example: a rule might be defined as:

    rule test:
        input: globspand('/my/path/in/{sample}/{batch}/{filename}.bam')
        output: '/my/path/out/{sample}/{batch}/result.tsv'

    In this case, the rule will run for exactly one sample and batch (the ones
    defined in the output) but will use all the bam files that complete the
    pattern in the input for that sample and batch.

    Alternatively, one can also apply static wildcards for expansion. For example:

    rule test_all:
        input: globspand('/my/path/in/{sample}/{batch}/{filename}.bam', sample=all_samples, batch=all_batches)
        output: '/my/path/out/all.tsv'

    Or they can even be used in combination (but must not define the same wildcards).
    This might be useful if you want to filter input files, like in this example:

    rule test_all:
        input: globspand('/my/path/in/{qc_type}/{sample}/{batch}/{filename}.bam', qc_types=['pass','fail'])
        output: '/my/path/out/{sample}/{batch}/result.tsv'
    """

    # All keys in the original path
    all_keys = [i[1] for i in Formatter().parse(path)  if i[1] is not None]
    # Keys provided statically during the invocation of globspand
    keys_in_globspand = kwargs.keys()
    def globspand_inner(wildcards):
        # Keys inferred from the output
        keys_in_wildcards = wildcards.keys()
        # Sanity check that there isn't an overlap:
        assert(len(set(keys_in_globspand).intersection(keys_in_wildcards))==0)
        # Keys that need to be found using glob_wildcard
        keys_glob = {k for k in all_keys if k not in keys_in_wildcards and k not in keys_in_globspand}
        # Expand the path using the wildcards from the output and the static wildcards
        # in globspand, while keeping a wildcard notation for those that aren't given
        # (that's wha the {{{k}}} does)
        wildcard_expand_dict = {k:(wildcards[k] if k in keys_in_wildcards
                                   else kwargs[k] if k in keys_in_globspand
                                   else '{{{k}}}'.format(k=k)) for k in all_keys}
        wildcard_expanded = expand(path, **wildcard_expand_dict)
        # Since we also allow lists in the static wildcards, this
        # might lead to multiple paths over which we iterate here
        for expanded_path in wildcard_expanded:
            # Find glob_wildcards for the expanded path and further expand it
            glob = glob_wildcards(expanded_path)
            glob_expand_dict = {k:getattr(glob,k) for k in keys_glob}
            globspanded_paths = expand(expanded_path, **glob_expand_dict)
            for globspanded_path in globspanded_paths:
                yield globspanded_path
    return globspand_inner


def makelist(x):
    """
    Utility function that returns a list [x] if x is a scalar and list(x)
    if x is iterable.
    """
    if isinstance(x, str):
        # If it's a string, make it a list
        return [x]
    try:
        iter(x)
    except TypeError:
        # If it's not iterable make it a list
        return [x]
    # Make sure iterable is a list
    return list(x)


def glob_output_from_input(in_pattern, out_pattern, combiner=zip, **wildcards_dict):
    """
    Another untility function that creates a list of input files based
    on the wildcards in an input file pattern and an output file pattern.
    The wildcards given in the wildcards_dict will be used to expand the patterns
    and all remaining wildcards in the output pattern will be used to glob_wildcard
    the input pattern to then further expand the output pattern.

    For example, if you had files

    /my/path/in/sample<n>/batch<m>/<l>.fq

    And a rule like

    rule test:
        input: globspand('/my/path/in/sample{n}/batch{m}/{l}.fq')
        output: '/my/path/out/sample{n}/{m}.bam'

    You could then run this rule for all test:

    rule all_test:
        input: glob_output_from_input('/my/path/in/sample{n}/batch{m}/{l}.fq', '/my/path/out/sample{n}/{m}.bam', n=[1,2])

    It would run all_test for samples 1 and 2, and find all the applicable batches n
    for which it can run test, and run all of those.
    """
    def glob_output_from_input_inner(_):
        in_keys = [i[1] for i in Formatter().parse(in_pattern)  if i[1] is not None]
        out_keys = [i[1] for i in Formatter().parse(out_pattern)  if i[1] is not None]

        unresolved_wildcards_output = [k for k in out_keys if k not in wildcards_dict.keys()]
        unresolved_wildcards_input = [k for k in in_keys if k not in wildcards_dict.keys()]
        assert(len(unresolved_wildcards_output)>0)
        wildcards_keys = [k for k in wildcards_dict.keys()]
        wildcards_values = [makelist(wildcards_dict[k]) for k in wildcards_keys]
        wildcards_values = combiner(*wildcards_values)
        for combination in wildcards_values:
            combination_dict = {wildcards_keys[i]:combination[i] for i in range(len(combination))}
            dummy_dict_in = {k:'{{{k}}}'.format(k=k) for k in unresolved_wildcards_input}
            dummy_dict_out = {k:'{{{k}}}'.format(k=k) for k in unresolved_wildcards_output}
            in_expanded = in_pattern.format(**combination_dict, **dummy_dict_in)
            out_expanded = out_pattern.format(**combination_dict, **dummy_dict_out)
            glob = glob_wildcards(in_expanded)
            glob_out_vals = []
            for key in unresolved_wildcards_output:
                glob_out_vals.append(getattr(glob, key))
            out_combinations_from_glob = list(set(zip(*glob_out_vals)))
            for out_combination in out_combinations_from_glob:
                out_combination_dict = {unresolved_wildcards_output[i]:out_combination[i] for i in range(len(out_combination))}
                out_formatted = out_expanded.format(**out_combination_dict)
                yield out_formatted
    return glob_output_from_input_inner

'''
##############################################################################
# Splitting fast5 files into batches for parallel processing computation
##############################################################################

The input directory raw/samplename/guppy/ will be searched for fast5 files.
The list of fast5 files per sample is split into evenly sized batches.

The directory raw/samplename/batched will contain a subdirectory for each 
batch.

CAVEAT: The fast5 files will be **symlinked** from the batch directory to the
basecalled directory. Do NOT delete the fast5 files from the basecalled
directory!

Rule will not be performed if "batched" directory exists. Delete the "batched"
directory, if you want to redo batch splitting.

This rule will be performed locally, as it is only creating symlinks and is
not computationally expensive.
'''


include: 'rules/fastq.rules'
include: 'rules/guppy.rules'
include: 'rules/mapping.rules'
include: 'rules/medaka.rules'
include: 'rules/megalodon.rules'
include: 'rules/nanopolish.rules'
include: 'rules/pycoqc.rules'
include: 'rules/sniffles.rules'
include: 'rules/tombo.rules'
include: 'rules/edgecase.rules'
include: 'rules/nanocompore.rules'
include: 'rules/albacore.rules'

# include: 'rules.d/custom.rules'

def split_batches_from_file_list(all_files, outdir):
    # Create "batched" directory if it doesn't exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    i = 0
    b = 0
    while i < len(all_files):
        batchdir = os.path.join(outdir, '%d' % b)
        os.mkdir(batchdir)
        for _ in range(per_batch):
            if i == len(all_files):
                break
            link = os.path.join(batchdir, '%s' % os.path.basename(all_files[i]))
            os.symlink(all_files[i], link)
            i += 1
        b += 1


def split_batches_input(wildcards):
    file_path = os.path.join(basedir, 'raw', wildcards['sample'], 'guppy', '{fname}.fast5')
    return expand(file_path, fname=glob_wildcards(file_path).fname)


checkpoint split_batches:
    input: split_batches_input
    output: os.path.join(basedir, 'raw', '{sample}', 'batched', 'done')
    params:
          outdir=os.path.join(basedir, 'raw', '{sample}', 'batched')
    run:
        split_batches_from_file_list(input, params.outdir)
        with open(output[0], 'w') as fp:
            pass


def split_batches_from_fastq_input(wildcards):
    file_path = os.path.join(basedir, 'fastq', wildcards['sample'], 'guppy', '{fname}.%s' % fastq_ending)
    return expand(file_path, fname=glob_wildcards(file_path.fname))


checkpoint split_batches_from_fastq:
    input: split_batches_from_fastq_input
    output: os.path.join(basedir, 'fastq', '{sample}', 'batched', 'done')
    params:
          outdir=os.path.join(basedir, 'fastq', '{sample}', 'batched')
    run:
        split_batches_from_file_list(input, params.outdir)
        with open(output[0], 'w') as fp:
            pass

checkpoint all_split_batches:
    input: expand(rules.split_batches.output, sample=unique_samples)

checkpoint all_split_batches_from_fastq:
    input: expand(rules.split_batches_from_fastq.output, sample=unique_samples)



'''
##############################################################################
# Reports methylation histograms - only really useful for benchmarking
# datasets where you know the methylation rate (like fully methylated or
# fully unmethylated datasets)
##############################################################################
'''


def report_methylation_input(wildcards):
    return expand(rules.metcall.output, batch=samplebatches[wildcards.sample],
                  sample=wildcards.sample, mtype=wildcards.mtype)


rule report_methylation:
    input: report_methylation_input
    output: os.path.join(basedir, 'report/{sample}_{mtype}.pdf')
    params:
          jobname='report_{sample}_{mtype}',
          runtime='0:10',
          memusage='4000',
          slots='1',
          misc=''
    shell: '{python} report_met.py ' + os.path.join(basedir, 'met') + \
         ' {output} {wildcards.sample} {wildcards.mtype}'

rule all_report_methylation:
    input: expand(rules.report_methylation.output, mtype=mettypes, sample=unique_samples)


localrules: prepare_mergebams, split_batches, split_batches_from_fastq
