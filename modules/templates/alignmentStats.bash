#!/usr/bin/env bash

set -euo pipefail
${params.alignmentStatsCommand} ${alignmentsSam} > ${sample}.alignmentStats.txt
