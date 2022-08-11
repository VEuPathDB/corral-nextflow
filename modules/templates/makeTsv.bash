#!/usr/bin/env bash

makeTsv.pl . .taxa.tsv ${params.summaryColumn} ${params.summaryFormat} > ${params.summaryColumn}.${params.summaryFormat}.tsv
