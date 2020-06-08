#!/bin/bash

# Setup some safe shell options
set -eu -o pipefail

# No arguments
if [ $# -lt 1 ]
then
	printf "\nPlease enter a case name from the examples...exiting\n\n"
	exit

# More than one argument
elif [ $# -gt 1 ]
then
	printf "\nToo many arguments...exiting\n\n"
	exit
fi

# Get case name
case=${1%/}

# Check if it is a real case
if [ ! -d $case ]; then
	printf "\nCase name is not a real example case...exiting\n\n"
	exit
fi

# Otherwise this a real case and ready to go
printf "\nRunning the $case example!\n\n"

# Clean the directory first
rm -rf $case/LIFE $case/Results

# Copy the params.h file into source folder
cp $case/params.h ../inc/params.h

# Build LIFE
(cd .. && make clean && make)

# Copy LIFE to case folder
cp ../LIFE $case/.

# Run LIFE
(cd $case && ./LIFE)

# Print finish
printf "\n\nFinished running the $case example!\n"

