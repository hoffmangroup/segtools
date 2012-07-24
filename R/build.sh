#!/usr/bin/env bash

## build.sh: build R package

## $Revision$
## Copyright 2012 Michael M. Hoffman <mmh1@uw.edu>

set -o nounset -o pipefail -o errexit

if [[ $# != 0 || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo usage: "$0"
    exit 2
fi

R CMD build segtools
