# Regression tests

Small, fast tests that check whether a code change alters SNEEP's output. They build
`differentialBindingAffinity_multipleSNPs` from `src/`, run it on a tiny synthetic data set and
compare all output files with stored reference outputs. A run takes a few seconds and needs
no data from Zenodo.

The tests check that the output **stays the same**, not that it is correct: the reference
outputs reflect the code at the time they were created, including known issues (see below).
Correctness is covered by the benchmarks in the paper and `runTests.sh` on real data.

## Run

```bash
bash tests/runRegressionTests.sh            # compare with tests/expected/<platform>/
CXX=clang++ bash tests/runRegressionTests.sh   # other compiler (default g++)
```

Requirements: a C++11 compiler, `bedtools`, `python3`. Without OpenMP a serial version is built
(the tests use `-n 1` anyway). If `python` is missing, the script provides a small wrapper (only
needed for SNEEP versions before 2026-10-07, whose Python helpers used `#!/usr/bin/env python`).

Output: `PASSED`, or the files that differ with the differing lines. All outputs and logs of the
last run are in `tests/out/` (not under version control).

## Test cases

| Case | Command (besides `-n 1 -p 0.5 -c 0.01 -b … -x …`) |
|---|---|
| `maf` | REMs (`-r`, `-g`), 3 background rounds matched by MAF (`-j 3 -k … -l 2`) |
| `gc` | as `maf`, matched by MAF and GC content (`-s true`) |
| `noREM` | no REMs, 2 background rounds (`-j 2 -l 5`) |
| `given` | random SNPs given (`-j 2 -i data/presampled`); results must equal `noREM`, and `data/presampled` must stay unchanged |

## Comparison (`compareOutputs.py`)

Every file of the run must exist in the reference and vice versa. Before comparing:
- lines are sorted (their order depends on threads and `std::sort`),
- `TF_count.txt` is transposed and sorted by TF (its column order follows the order of the motif
  files in the directory, which depends on the file system),
- the date line of `info.txt` is dropped,
- numbers may differ by a relative tolerance of 1e-5 (last printed digit, e.g. other math libraries).

## Reference outputs per platform

`tests/expected/<uname -s>/` (`Darwin`, `Linux`). The random numbers of the background sampling
(`std::uniform_int_distribution`) differ between libc++ (macOS) and libstdc++ (Linux), so the
sampled SNVs differ between platforms. The first run on a new platform creates its reference with
`--update` (only from a code version whose output is trusted).

## Intended changes of the output

If a change is supposed to alter the output (e.g. a bug fix), run the tests, check that only the
expected lines differ, then store the new reference and commit it together with the code change:

```bash
bash tests/runRegressionTests.sh --update
```

## Test data (`data/`)

Created by `makeTestData.py` (deterministic, fixed seed); only needed to recreate or extend the data:
`python3 tests/makeTestData.py . tests/data` (from `SNEEP-software/`), then `--update`.

- `genome.fa`: synthetic chr1 (20 kb) and chr2 (15 kb); planted instances of CTCF(MA0139.1), GATA1,
  FOXA1 and SP1 with SNVs causing loss and gain of binding, also on the reverse strand; the GATA1
  instances lie in a soft-masked (lower case) stretch; N blocks on chr2
- `motifs.transfac`, `scales.txt`: the 4 motifs and scales from the real files
- `snvs.bed`: 26 input SNVs, plus a duplicate line and an indel
- `interactions.txt`, `geneNames.txt`: 6 REMs of 6 genes, two of them overlapping
- `presampled/`: `randomSNPs_0.txt` and `randomSNPs_1.txt`, a copy of the random SNPs of case `noREM`
  (input for `-i`; not created by `makeTestData.py`, copy them again if the data set changes)
- `dbSNPs_sorted.txt`: about 3,500 SNVs in the format of `getSNPInfo` (incl. GC content, column 9),
  with MAF -1, MAF 0 and multi-allelic entries

Expected behaviour visible in the reference: CTCF, FOXA1 and SP1 sites are significant; GATA1 is
just not significant (p ≈ 0.011); the SNV between the N blocks triggers the GC warning; in the `gc`
case some GC bins fall back to the nearest bin (warnings in `tests/out/logs/gc.txt`).

## Known issues contained in the reference (2026-10-07)

- Very strong hits get p-value 0 (`1 - (1 - e^(-x/b))^n` underflows below ~1e-16), e.g. the CTCF SNVs.

Fixing them will change the output; update the reference then.
