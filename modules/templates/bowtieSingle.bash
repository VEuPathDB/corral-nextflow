#!/usr/bin/env bash

set -euo pipefail
grep -c '^@' ${readsFastq} > numReads.txt
${params.bowtie2Command} \
  -x ${params.refdb} \
  -U ${readsFastq} \
  -S alignmentsSingle.sam 
