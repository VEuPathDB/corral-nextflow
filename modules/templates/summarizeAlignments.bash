#!/usr/bin/env bash

${params.summarizeAlignmentsCommand} \
  --input ${alignmentsSam} \
  --refdb-marker-to-taxon-path ${params.markerToTaxonPath} \
  --refdb-format eukprot \
  --output-type taxon_all \
  --num-reads \$(cat ${numReadsPath}) \
  --output ${sample}.taxa.tsv 
