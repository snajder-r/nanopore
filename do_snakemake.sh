#!/bin/bash

cluster_script="bsub {params.misc} -n {params.slots} -W {params.runtime} -M {params.memusage}  -R \"rusage[mem={params.memusage}]\" -o {basedir}/logs/{params.jobname}.log -e {basedir}/logs/{params.jobname}.err"

snakemake --cluster "$cluster_script" --jobs 128 --latency-wait 120 $@

