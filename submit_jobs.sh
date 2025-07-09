#!/bin/bash

rm -rf "results"

mkdir -p results
bsub <./runjobs.lsf 


echo "DONE Submitting Jobs."
