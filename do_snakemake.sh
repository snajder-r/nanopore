#!/bin/bash

if [[ -z "${MEDUSA_CLUSTER_SCRIPT}" ]]; then
  echo "Please set environment variable MEDUSA_CLUSTER_SCRIPT." >&2;
  exit 1
else
  cluster_script="${MEDUSA_CLUSTER_SCRIPT}"
fi

if [[ -z "${HDF5_PLUGIN_PATH}" ]]; then
  echo "Please set environment variable HDF5_PLUGIN_PATH to a path that contains the vbz plugin. See https://github.com/nanoporetech/vbz_compression." >&2;
  exit 1
fi

snakemake --cluster "$cluster_script" --jobs 128 --latency-wait 120 $@

