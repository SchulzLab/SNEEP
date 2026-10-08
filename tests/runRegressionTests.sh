#!/usr/bin/env bash
# Regression tests for SNEEP: builds the tools from src/, runs them on the small data set in tests/data/
# and compares the outputs with tests/expected/<platform>/ (see tests/README.md).
#
# Usage: bash tests/runRegressionTests.sh [--update]
#	--update	store the current outputs as the new expected outputs (only after checking the differences!)
# Environment: CXX (default g++), needs bedtools and python3 in PATH

set -u
TESTS=$(cd "$(dirname "$0")" && pwd)
SRC=$TESTS/../src
WORK=$TESTS/out
PLATFORM=$(uname -s) # one reference per platform: random numbers (std::uniform_int_distribution) differ between standard libraries
EXPECTED=$TESTS/expected/$PLATFORM
CXX=${CXX:-g++}

rm -rf "${WORK:?}"
mkdir -p "$WORK/bin" "$WORK/results" "$WORK/logs"

# --- build (OpenMP if available, otherwise serial with a stub omp.h)
if echo 'int main(){}' | $CXX -x c++ -fopenmp - -o "$WORK/bin/omptest" 2>/dev/null; then
	OMP="-fopenmp"
else
	mkdir -p "$WORK/stub"
	printf '#pragma once\ninline int omp_get_thread_num(){return 0;}\ninline int omp_get_num_threads(){return 1;}\n' > "$WORK/stub/omp.h"
	OMP="-I$WORK/stub"
	echo "note: $CXX has no OpenMP, building serial versions"
fi
for tool in differentialBindingAffinity_multipleSNPs; do
	if ! $CXX -std=c++11 -O2 $OMP "$SRC/$tool.cpp" -o "$WORK/bin/$tool" 2> "$WORK/logs/build_$tool.txt"; then
		echo "build of $tool failed, see $WORK/logs/build_$tool.txt"
		exit 1
	fi
done
# fingerprint of the random numbers of this compiler/standard library (see rngFingerprint.cpp)
$CXX -std=c++11 "$TESTS/rngFingerprint.cpp" -o "$WORK/bin/rngFingerprint" 2> "$WORK/logs/build_rngFingerprint.txt" || { echo "build of rngFingerprint failed"; exit 1; }
FINGERPRINT=$("$WORK/bin/rngFingerprint" | cksum | awk '{print $1}')
# the python helpers are called via PATH; their shebang is "python", which might only exist as python3
if ! python -c "" > /dev/null 2>&1; then # macOS has a stub "python" that only asks to install the developer tools
	mkdir -p "$WORK/python"
	printf '#!/bin/sh\nexec python3 "$@"\n' > "$WORK/python/python" # no symlink: macOS python3 is a launcher that checks its own name
	chmod +x "$WORK/python/python"
	export PATH="$WORK/python:$PATH"
fi
export PATH="$WORK/bin:$SRC:$PATH"

# --- run (relative paths, so the paths written to info.txt do not depend on the location of the repository)
cd "$TESTS"
R=out/results
run(){ # name, command...
	local name=$1
	shift
	echo "run $name"
	"$@" > "out/logs/$name.txt" 2>&1 || echo "  exit code $? (see tests/out/logs/$name.txt)"
}
COMMON="-n 1 -p 0.5 -c 0.01 -b ../necessaryInputFiles/frequency.txt -x ../necessaryInputFiles/transition_matrix.txt"
INPUT="data/motifs.transfac data/snvs.bed data/genome.fa data/scales.txt"
DBSNP=data/dbSNPs_sorted.txt # mini version of the file written by getSNPInfo (incl. GC content)

# 1. REMs, background sampling matched by MAF
run maf differentialBindingAffinity_multipleSNPs -o $R/maf/ $COMMON -j 3 -k $DBSNP -l 2 -r data/interactions.txt -g data/geneNames.txt $INPUT
# 2. same, matched by MAF and GC content
run gc differentialBindingAffinity_multipleSNPs -o $R/gc/ $COMMON -j 3 -k $DBSNP -l 2 -s true -r data/interactions.txt -g data/geneNames.txt $INPUT
# 3. without REMs
run noREM differentialBindingAffinity_multipleSNPs -o $R/noREM/ $COMMON -j 2 -k $DBSNP -l 5 $INPUT
# 4. random SNPs given (-i), read-only use of data/presampled (copy of the random SNPs of case noREM)
PRESAMPLED_BEFORE=$(cat data/presampled/* | cksum)
run given differentialBindingAffinity_multipleSNPs -o $R/given/ $COMMON -j 2 -i data/presampled $INPUT
if [ "$(ls data/presampled | wc -l | tr -d ' ')" != "2" ] || [ "$(cat data/presampled/* | cksum)" != "$PRESAMPLED_BEFORE" ]; then
	echo "FAILED: SNEEP changed the -i directory data/presampled"
	exit 1
fi

# --- compare or update
if [ "${1:-}" == "--update" ]; then
	rm -rf "${EXPECTED:?}"
	mkdir -p "$EXPECTED"
	cp -R "$WORK/results/." "$EXPECTED/"
	echo "$FINGERPRINT" > "$EXPECTED/rng_fingerprint.txt"
	echo "expected outputs updated: $EXPECTED"
	exit 0
fi
if [ ! -d "$EXPECTED" ]; then
	echo "no expected outputs for $PLATFORM yet; create them with: bash tests/runRegressionTests.sh --update"
	exit 1
fi
# sampling-dependent outputs are only comparable with the same random numbers (same compiler/standard library)
MODE=""
if [ ! -f "$EXPECTED/rng_fingerprint.txt" ]; then
	echo "note: no random number fingerprint stored with the reference, comparing all outputs (store it with --update)"
elif [ "$(cat "$EXPECTED/rng_fingerprint.txt")" != "$FINGERPRINT" ]; then
	echo "note: this compiler/standard library draws other random numbers than the one of the reference (std::uniform_int_distribution),"
	echo "      so the sampled random SNPs differ; comparing only the outputs that do not depend on them (input SNPs, case given)"
	MODE="--without-sampling"
fi
python3 "$TESTS/compareOutputs.py" $MODE "$EXPECTED" "$WORK/results" && echo "PASSED" || { echo "FAILED"; exit 1; }
