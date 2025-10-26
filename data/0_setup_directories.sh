#!/usr/bin/env bash

#====================================================#
# * Purpose: To create a file structure that makes sense
# for the data we will be using
#
# * Requires: "the-milk-man" conda env I made for this; I
# will make the yaml availble on the GitHub

# i wanted to create more robust bash scripts. there are a 
# lot of fail conditions that fuck them up. credit to this 
# gist for explaining some of the best conditions [s]. i 
# could of also added '-x' to debug print to terminal.
#
# [s]: https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -ueo pipefail

#===== Setup
printf "(1) Checking if conda is activate...\n"
if [[ $CONDA_DEFAULT_ENV == "the-milk-man" ]]; then
	printf "\t--> The conda env you need is already activate, good...\n\n"
elif [[ $CONDA_DEFAULT_ENV != "the-milk-man" ]]; then
	printf "\t--> You need to activate 'the-milk-man'.\n\n"
	exit 1
else
	printf "\t--> Something is wrong...\n\n"
	exit 1
fi

printf "(2) Creating storage directories for the data...\n"
if [ -d controls ]; then
	printf "\t--> The control directory exists, good.\n"
else
	printf "\t--> Creating the control directory!\n"
	mkdir controls
fi

if [ -d experiments/peak ]; then
	printf "\t--> The peak experiments directory exists, good.\n"
else
	printf "\t--> Creating the peak experiments directory!\n"
	mkdir -p experiments/peak
fi

if [ -d experiments/late ]; then
	printf "\t--> The late experiments directory exists, good.\n\n"
else
	printf "\t--> Creating the late experiments directory!\n\n"
	mkdir -p experiments/late
fi

#===== Arguments
LIST_EXP_PEAK="$(pwd)/experiments-peak.txt"
LIST_EXP_LATE="$(pwd)/experiments-late.txt"
LIST_CTRLS="$(pwd)/controls.txt"

#===== Pipline

# manage controls first
for i in $(cat controls.txt); do prefetch --progress $i; done
for i in $(cat controls.txt); do fasterq-dump $i --split-files -p -e 8 -O controls/; done
pigz -p 8 controls/*

# then peak files
for i in $(cat experiments-peak.txt); do prefetch --progress $i; done
for i in $(cat experiments-peak.txt); do fasterq-dump $i --split-files -p -e 8 -O experiments/peak/; done
pigz -p 8 experiments/peak/*

# finally, late files
for i in $(cat experiments-late.txt); do prefetch --progress $i; done
for i in $(cat experiments-late.txt); do fasterq-dump $i --split-files -p -e 8 -O experiments/late/; done
pigz -p 8 experiments/late/*

exit 0
