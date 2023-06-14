#!/usr/bin/env bash

set -euo pipefail

if [ -s $HOME/.ncbi/user-settings.mkfg ]
then
    fasterq-dump --split-3 ${id}
else
    mkdir -p $HOME/.ncbi/
    fasterq-dump --split-3 ${id}
fi
