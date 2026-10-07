#!/usr/bin/env python3
"""
Creates the small synthetic data set for the regression tests in tests/data/.

Not needed to run the tests (the data is checked in); only to recreate or extend
the data set. Deterministic (fixed seed), so rerunning it gives the same files.

Usage: python3 makeTestData.py <SNEEP-software dir> <output dir>

Creates:
	genome.fa			two synthetic chromosomes (chr1 20 kb, chr2 15 kb) with
					planted motif instances, a soft-masked (lower case) stretch
					and N blocks
	motifs.transfac			4 motifs taken from examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt
	scales.txt			their scales from necessaryInputFiles/estimatedScalesPerMotif_1.9.txt
	snvs.bed			input SNVs (SNEEP bed-like format), incl. a duplicate and an indel
	interactions.txt		REMs linked to genes (-r), 12 columns
	geneNames.txt			ensemblID,geneName (-g)
	dbSNPs_sorted.txt		mini dbSNP file in the format written by getSNPInfo (-k), incl. GC content
"""

import os
import random
import sys

MOTIFS = ["CTCF(MA0139.1)", "GATA1", "FOXA1", "SP1"]
CHROM_LENGTHS = {"chr1": 20000, "chr2": 15000}
GC_FLANK = 30  # as in getSNPInfo.cpp
BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def read_motifs(transfac_file):
	"""returns {ID: (transfac block, count matrix as list of [A, C, G, T])}"""
	motifs = {}
	block = []
	for line in open(transfac_file):
		block.append(line)
		if line.startswith("//"):
			name = [l.split()[1] for l in block if l.startswith("ID")][0]
			counts = []
			in_matrix = False
			for l in block:
				if l.startswith("PO"):
					in_matrix = True
				elif in_matrix and l[:2].isdigit():
					counts.append([float(x) for x in l.split()[1:5]])
				elif in_matrix:
					in_matrix = False
			motifs[name] = ("".join(block), counts)
			block = []
	return motifs


def reverse_complement(seq):
	return "".join(COMPLEMENT[b] for b in reversed(seq))


def main():
	if len(sys.argv) < 3:
		print("Usage: python3 makeTestData.py <SNEEP-software dir> <output dir>")
		sys.exit(1)
	sneep_dir, out_dir = sys.argv[1], sys.argv[2]
	os.makedirs(out_dir, exist_ok=True)
	rng = random.Random(1)

	# motifs and scales
	motifs = read_motifs(os.path.join(sneep_dir, "examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt"))
	with open(os.path.join(out_dir, "motifs.transfac"), "w") as f:
		for m in MOTIFS:
			f.write(motifs[m][0])
	with open(os.path.join(sneep_dir, "necessaryInputFiles/estimatedScalesPerMotif_1.9.txt")) as f, \
			open(os.path.join(out_dir, "scales.txt"), "w") as out:
		lines = f.readlines()
		out.write(lines[0])
		for line in lines[1:]:
			if line.split("\t")[0] in MOTIFS:
				out.write(line)

	# random genome (AT rich like the human genome)
	genome = {c: [rng.choices(BASES, weights=[29, 21, 21, 29])[0] for _ in range(n)] for c, n in CHROM_LENGTHS.items()}

	snvs = []  # (chr, pos, var1, var2, rsID, MAF)
	maf_values = ["0.012", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "1e-07", "-1"]

	# planted motif instances on chr1: per motif 3 instances with a SNV at an informative position
	# 1: loss of binding (var1 = consensus base), 2: gain of binding (genome carries the worst base),
	# 3: reverse strand, loss at the second most informative position
	# the instances of the second motif lie in a soft-masked (lower case) stretch
	start = 1000
	for k, m in enumerate(MOTIFS):
		counts = motifs[m][1]
		consensus = "".join(BASES[row.index(max(row))] for row in counts)
		information = sorted(range(len(counts)), key=lambda i: -max(counts[i]) / sum(counts[i]))
		for inst in range(3):
			seq = list(consensus)
			pos_in_motif = information[0] if inst < 2 else information[1]
			best = consensus[pos_in_motif]
			worst = BASES[counts[pos_in_motif].index(min(counts[pos_in_motif]))]
			if inst == 1:
				seq[pos_in_motif] = worst
			site = "".join(seq)
			if inst == 2:
				site = reverse_complement(site)
				pos_in_site = len(site) - 1 - pos_in_motif
				best, worst = COMPLEMENT[best], COMPLEMENT[worst]
			else:
				pos_in_site = pos_in_motif
			s = start + inst * 400
			genome["chr1"][s:s + len(site)] = list(site)
			pos = s + pos_in_site
			if inst == 1:
				var1, var2 = worst, best  # gain
			else:
				var1, var2 = best, worst  # loss
			snvs.append(("chr1", pos, var1, var2, "rsMotif%d_%d" % (k, inst), maf_values[(3 * k + inst) % len(maf_values)]))
		start += 1500
	for i in range(2500, 3800):  # soft-masked stretch (instances of the second motif)
		genome["chr1"][i] = genome["chr1"][i].lower()

	# N blocks on chr2 around one SNV (GC window nearly only N) and a block next to another SNV
	for i in list(range(7000, 7070)) + list(range(7071, 7141)):
		genome["chr2"][i] = "N"
	for i in range(9000, 9040):
		genome["chr2"][i] = "N"
	for chrom, pos in (("chr2", 7070), ("chr2", 9045)):
		ref = genome[chrom][pos].upper()
		snvs.append((chrom, pos, ref, rng.choice([b for b in BASES if b != ref]), "rsN_%d" % pos, "0.3"))

	# random SNVs
	for i in range(12):
		chrom = "chr1" if i < 6 else "chr2"
		while True:
			pos = rng.randint(200, CHROM_LENGTHS[chrom] - 200)
			if chrom == "chr2" and (6900 < pos < 7250 or 8900 < pos < 9150):
				continue
			if chrom == "chr1" and 900 < pos < 7200:
				continue
			break
		ref = genome[chrom][pos].upper()
		snvs.append((chrom, pos, ref, rng.choice([b for b in BASES if b != ref]), "rsRandom%d" % i, maf_values[i % len(maf_values)]))

	with open(os.path.join(out_dir, "genome.fa"), "w") as f:
		for chrom, seq in genome.items():
			f.write(">" + chrom + "\n")
			seq = "".join(seq)
			for i in range(0, len(seq), 60):
				f.write(seq[i:i + 60] + "\n")

	with open(os.path.join(out_dir, "snvs.bed"), "w") as f:
		for chrom, pos, var1, var2, rs, maf in snvs:
			f.write("\t".join([chrom, str(pos), str(pos + 1), var1, var2, rs, maf]) + "\n")
		chrom, pos, var1, var2, rs, maf = snvs[0]  # duplicate line -> removed by SNEEP
		f.write("\t".join([chrom, str(pos), str(pos + 1), var1, var2, rs, maf]) + "\n")
		f.write("\t".join(["chr1", "15000", "15001", "A", "AT", "rsIndel", "0.2"]) + "\n")  # indel -> not considered

	# REMs: two overlapping REMs (two genes) around the instances of the first motif, one per other motif, one on chr2
	rems = [("chr1", 900, 2000, 1), ("chr1", 1300, 2600, 2), ("chr1", 2400, 3900, 3),
		("chr1", 3900, 5400, 4), ("chr1", 5400, 6900, 5), ("chr2", 200, 15000, 6)]
	with open(os.path.join(out_dir, "interactions.txt"), "w") as f:
		for chrom, s, e, g in rems:
			f.write("\t".join([chrom, str(s), str(e), "ENSG0000000000%d.1" % g, "REM%07d" % g,
				"%.3f" % rng.uniform(-1, 1), "%.2e" % rng.uniform(1e-6, 1e-2), ".", ".", ".", "-", "-"]) + "\n")
	with open(os.path.join(out_dir, "geneNames.txt"), "w") as f:
		for g in range(1, 7):
			f.write("ENSG0000000000%d.1,GENE%d\n" % (g, g))

	# mini dbSNP file in the format of getSNPInfo (dbSNPs_sorted.txt): MAF chr start end ref alt rsID MAF GC,
	# sorted by MAF; MAF uniform in [0, 0.5], some with MAF 0 (as in the real dbSNP), some without MAF (-1),
	# some multi-allelic; GC as in getSNPInfo: (C + G) / (A + C + G + T) in +- 30 bp, N excluded, 6 digits
	dbsnp = []
	for chrom, n in CHROM_LENGTHS.items():
		seq = "".join(genome[chrom]).upper()
		for pos in sorted(rng.sample(range(150, n - 150), 2000 if chrom == "chr1" else 1500)):
			ref = seq[pos]
			if ref == "N":
				continue
			alts = [b for b in BASES if b != ref]
			alt = ",".join(rng.sample(alts, 2)) if rng.random() < 0.05 else rng.choice(alts)
			r = rng.random()
			if r < 0.1:
				maf = "-1"
			elif r < 0.13:
				maf = "0"
			else:
				maf = "%g" % round(rng.uniform(0.0005, 0.5), 4)
			window = seq[max(0, pos - GC_FLANK):pos + GC_FLANK + 1]
			acgt = sum(window.count(b) for b in BASES)
			gc = "%.6g" % ((window.count("C") + window.count("G")) / acgt) if acgt > 0 else "-1"
			dbsnp.append((float(maf), "\t".join([maf, chrom, str(pos), str(pos + 1), ref, alt, "rs%s%d" % (chrom[-1], pos), maf, gc])))
	dbsnp.sort(key=lambda x: x[0])
	with open(os.path.join(out_dir, "dbSNPs_sorted.txt"), "w") as f:
		for _, line in dbsnp:
			f.write(line + "\n")


if __name__ == "__main__":
	main()
