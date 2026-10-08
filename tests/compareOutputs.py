#!/usr/bin/env python3
"""
Compares the outputs of a regression test run with the expected outputs.

Usage: python3 compareOutputs.py <expected dir> <actual dir>

Every file below <expected dir> must exist below <actual dir> with the same content, and vice versa.
Normalisation before comparing (to be independent of thread scheduling and file system order):
	- lines are sorted (the order of SNVs/TFs in the outputs depends on threads and std::sort)
	- TF_count.txt is transposed and sorted by TF (column order = order of the motif files in the directory)
	- info.txt: the line with date and time is dropped
	- FASTA headers: the suffix "::chr:start-end" is dropped (added by newer bedtools versions only; SNEEP ignores it)
	- numbers may differ by a relative tolerance of 1e-5 (last printed digit, e.g. different math libraries)
Exit code 0 if all files match, 1 otherwise.
"""

import difflib
import os
import sys

REL_TOL = 1e-5
MAX_SHOWN = 10


def normalise(path):
	with open(path) as f:
		lines = f.read().splitlines()
	name = os.path.basename(path)
	lines = [l.split("::")[0] if l.startswith(">") else l for l in lines]
	if name == "info.txt":
		lines = [l for l in lines if not l.startswith("#\tdate and time")]
	if name == "TF_count.txt" and lines:
		rows = [l.split("\t") for l in lines]
		lines = ["\t".join(col) for col in zip(*rows)]
	return sorted(lines)


def same_token(a, b):
	if a == b:
		return True
	try:
		x, y = float(a), float(b)
	except ValueError:
		return False
	return abs(x - y) <= REL_TOL * max(abs(x), abs(y)) or abs(x - y) < 1e-300


def same_line(a, b):
	ta, tb = a.split("\t"), b.split("\t")
	return len(ta) == len(tb) and all(same_token(x, y) for x, y in zip(ta, tb))


def files_below(directory):
	result = set()
	for root, _, files in os.walk(directory):
		for f in files:
			result.add(os.path.relpath(os.path.join(root, f), directory))
	return result


def main():
	if len(sys.argv) < 3:
		print("Usage: python3 compareOutputs.py <expected dir> <actual dir>")
		sys.exit(2)
	expected_dir, actual_dir = sys.argv[1], sys.argv[2]
	expected, actual = files_below(expected_dir), files_below(actual_dir)
	failed = 0
	for f in sorted(expected - actual):
		print("MISSING  " + f)
		failed += 1
	for f in sorted(actual - expected):
		print("EXTRA    " + f)
		failed += 1
	for f in sorted(expected & actual):
		e, a = normalise(os.path.join(expected_dir, f)), normalise(os.path.join(actual_dir, f))
		if e == a:
			continue
		if len(e) == len(a) and all(same_line(x, y) for x, y in zip(e, a)):
			print("OK (numbers within tolerance)  " + f)
			continue
		failed += 1
		print("DIFFERS  " + f)
		diff = [l for l in difflib.unified_diff(e, a, "expected", "actual", lineterm="", n=0) if not l.startswith("@@")]
		for line in diff[:MAX_SHOWN + 2]:
			print("    " + line[:200])
		if len(diff) > MAX_SHOWN + 2:
			print("    ... %d more lines" % (len(diff) - MAX_SHOWN - 2))
	print("%d of %d files differ" % (failed, len(expected | actual)))
	sys.exit(1 if failed else 0)


if __name__ == "__main__":
	main()
