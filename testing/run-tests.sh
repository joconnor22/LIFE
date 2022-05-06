#!/bin/bash

# Setup some safe shell options
set -eu -o pipefail
shopt -s nullglob

if [ ! -d ../examples ]; then
    echo "../examples not found.  Please run run-tests.sh from the 'testing' directory." >&2
    exit 1
fi

# No arguments
if [ $# -eq 0 ]
then

	# Set testCases
	testCases=$(ls -d -- */)

# One argument
elif [ $# -eq 1 ]
then

	# Check if it is a real case
	if [ ! -d "$1" ]; then
		printf "\nCase name is not a real example case...exiting\n\n"
		exit
	fi

	# Set testCases
	testCases=$1

# More than one argument
elif [ $# -gt 1 ]
then
	printf "\nToo many arguments...exiting\n\n"
	exit
fi

# Run all the cases
for d in $testCases
do

	# Get case name
	case=${d%/}

	# Print header
	printf "\nRunning %s test!\n\n" "$case"

	# Check if there is reference data
	if [ ! -d "$case"/RefData ]; then
		printf "\nThere is no reference data for %s case...exiting\n\n" "$case"
		exit
	fi

	# Clean the directory first
	rm -rf "$case"/LIFE "$case"/Results "$case/$case.diff"

	# Copy params.h from RefData into src folder
	cp "$case"/RefData/params.h ../inc/params.h

	# Build LIFE
	(cd .. && make clean && make -j 8)

	# Copy case to directory
	cp ../LIFE "$case/."

	# Check if there a geometry.config file
	if [ -d "$case"/RefData/input ]; then
		cp -r "$case"/RefData/input "$case"/.
	fi

	# Run the case
	(cd "$case" && ./LIFE)

	# If TurekHron case then run again to test restart feature
	if [ "$case" == "TurekHron" ]; then
		(cd "$case" && ./LIFE)
	fi

	# Print finish
	printf "Finished running %s test!\n\n" "$case"
done

# Print header
printf "\n\n\nChecking against reference data...\n\n"

# Set colors
normal=$(tput sgr0)
red=$(tput setaf 1)
green=$(tput setaf 2)

# Number of fails
nFails=0

# Check the output against reference data
for d in $testCases
do

	# Get case name
	case=${d%/}

	# Check if results exist
	if [ -d "$case"/Results ]; then

		# Run diff and output to file
		if diff -r "$case"/Results "$case"/RefData/Results --exclude=Log.out &> "$case/$case.diff"; then
			printf "%s\n\n" "${green}PASS${normal} -> $case"
		else
			printf "%s\n\n" "${red}FAIL${normal} -> $case"
			nFails=$((nFails+1))
		fi
	else
		printf "%s\n\n" "${red}MISSING${normal} -> $case"
		nFails=$((nFails+1))
	fi
done

# Check if fails is more than zero times
if [ $nFails -eq 0 ]; then
	printf "\n%s\n\n" "${green}PASSED ALL TESTS!${normal}"
else
	printf "\n%s\n\n" "${red}FAILED SOME TESTS!${normal}"
fi
