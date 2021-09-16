#!/bin/bash

# Snakemake expects the submit script to print only the jobid, so we need this wrapper
bsub $@ | sed 's/.*<([0-9]+)>.*'
