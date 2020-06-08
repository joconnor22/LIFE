#!/bin/bash

# Setup some safe shell options
set -eu -o pipefail

# Find all example cases
for d in ../examples/*/
do

	# Get case path
	casePath=${d%/}

	# Get case name
	caseName=${casePath##*/}

	# Print header
	printf "\nStoring new $caseName RefData!\n\n"

	# Clean/create the case directory
	rm -rf $caseName
	mkdir $caseName

	# Copy params.h from examples into src folder
	cp $casePath/params.h ../inc/params.h

	# Modifiy the write out frequencies for testing
	sed -i "/const int nSteps/c\const int nSteps = 500;" ../inc/params.h
	sed -i "/const int tinfo/c\const int tinfo = nSteps / 50;" ../inc/params.h
	sed -i "/const int tVTK/c\const int tVTK = nSteps / 10;" ../inc/params.h
	sed -i "/const int tRestart/c\const int tRestart = nSteps / 5;" ../inc/params.h

	# Build LIFE
	(cd .. && make clean && make -j 8)

	# Create ref directory
	mkdir $caseName/RefData

	# Copy case to RefData
	cp ../LIFE $caseName/RefData/.
	cp ../inc/params.h $caseName/RefData/.

	# Check if there a geometry.config file
	if [ -d $casePath/input ]; then
		cp -r $casePath/input $caseName/RefData/.
	fi

	# Run the case
	(cd $caseName/RefData && ./LIFE)

	# If TurekHron case then run again to test restart feature
	if [ $caseName == "TurekHron" ]; then
		(cd $caseName/RefData && ./LIFE)
	fi

	# Print finish
	printf "Finished storing new $caseName RefData!\n\n"
done
