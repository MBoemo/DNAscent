#----------------------------------------------------------
# Copyright 2026 University of Cambridge
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


#--------------------------------------------------------------------------------------------------------------------------------------
def parseSeeBreaks(path):

	header = {}
	sections = {}
	current = None

	with open(path, "r") as fh:
		for raw in fh:
			line = raw.rstrip("\n")
			if not line:
				continue
			if line.startswith("#"):
				parts = line[1:].split(None, 1)
				key = parts[0]
				val = parts[1] if len(parts) > 1 else ""
				header[key] = val
				current = None
			elif line.startswith(">"):
				current = line[1:].rstrip(":")
				sections[current] = []
			elif current is not None:
				sections[current].append(line)

	return header, sections


#--------------------------------------------------------------------------------------------------------------------------------------
def asFloats(sections, name):
	return np.array([float(x) for x in sections.get(name, [])], dtype=float)


#--------------------------------------------------------------------------------------------------------------------------------------
def getFloat(header, key, default=float("nan")):
	try:
		return float(header[key])
	except (KeyError, ValueError):
		return default


#--------------------------------------------------------------------------------------------------------------------------------------
def plotEffect(header, sections, outPath):

	nullEffect = asFloats(sections, "NullEffectDistribution")
	effectBoot = asFloats(sections, "ObservedEffectBootstrap")

	if nullEffect.size == 0 or effectBoot.size == 0:
		sys.stderr.write("Skipping effect plot: NullEffectDistribution / ObservedEffectBootstrap not found.\n")
		return

	difMean = getFloat(header, "Difference")
	lowerBound = getFloat(header, "OneSided95LowerBound")
	pValue = getFloat(header, "OneSidedPValue")

	allVals = np.concatenate([nullEffect, effectBoot])
	bins = np.linspace(allVals.min(), allVals.max(), 60)

	fig, ax = plt.subplots(figsize=(7, 4.5))
	ax.hist(nullEffect, bins=bins, density=True, alpha=0.55, color="#9e9e9e", label="Null (random placement)")
	ax.hist(effectBoot, bins=bins, density=True, alpha=0.55, color="#1f77b4", label="Observed estimate (bootstrap)")

	if np.isfinite(difMean):
		ax.axvline(difMean, color="#1f77b4", linestyle="-", linewidth=1.5, label="Effect estimate")
	ax.axvline(0.0, color="#000000", linestyle=":", linewidth=1.2, label="No effect")
	if np.isfinite(lowerBound):
		ax.axvline(lowerBound, color="#d62728", linestyle="--", linewidth=1.5, label="One-sided 95% lower bound")

	# Shade just the null-density tail at or above the observed effect: this area is the p-value
	if np.isfinite(difMean):
		nullCounts, nullEdges = np.histogram(nullEffect, bins=bins, density=True)
		edgesTail = [difMean]
		heights = []
		for i in range(nullCounts.size):
			if nullEdges[i + 1] <= difMean:
				continue
			heights.append(nullCounts[i])
			edgesTail.append(nullEdges[i + 1])
		if heights:
			x = np.array(edgesTail)
			y = np.array(heights + [heights[-1]])
			ax.fill_between(x, y, step="post", color="#d62728", alpha=0.3, label="p-value (null tail)")

	title = "Excess analogue tracks at read ends"
	if np.isfinite(pValue):
		title += "   (one-sided p = {:.3g})".format(pValue)
	ax.set_title(title)
	ax.set_xlabel("Effect: observed - expected read-end fraction")
	ax.set_ylabel("Density")
	ax.legend(fontsize=8, loc="upper right", alpha=0.3)
	fig.tight_layout()
	fig.savefig(outPath, dpi=200)
	plt.close(fig)


#--------------------------------------------------------------------------------------------------------------------------------------
def plotRawFractions(header, sections, outPath):

	expected = asFloats(sections, "ExpectedReadEndFractions")
	observed = asFloats(sections, "ObservedReadEndFractions")

	if expected.size == 0 or observed.size == 0:
		sys.stderr.write("Skipping raw-fraction plot: ExpectedReadEndFractions / ObservedReadEndFractions not found.\n")
		return

	allVals = np.concatenate([expected, observed])
	bins = np.linspace(allVals.min(), allVals.max(), 60)

	fig, ax = plt.subplots(figsize=(7, 4.5))
	ax.hist(expected, bins=bins, density=True, alpha=0.55, color="#9e9e9e", label="Expected (null)")
	ax.hist(observed, bins=bins, density=True, alpha=0.55, color="#1f77b4", label="Observed")

	expMean = getFloat(header, "ExpectedReadEndFraction")
	obsMean = getFloat(header, "ObservedReadEndFraction")
	if np.isfinite(expMean):
		ax.axvline(expMean, color="#9e9e9e", linestyle="--", linewidth=1.5)
	if np.isfinite(obsMean):
		ax.axvline(obsMean, color="#1f77b4", linestyle="--", linewidth=1.5)

	ax.set_title("Raw read-end fractions (QC, pooled over tolerances)")
	ax.set_xlabel("Read-end fraction")
	ax.set_ylabel("Density")
	ax.legend(fontsize=8, loc="upper right", alpha=0.3)
	fig.tight_layout()
	fig.savefig(outPath, dpi=200)
	plt.close(fig)


#--------------------------------------------------------------------------------------------------------------------------------------
def main():

	parser = argparse.ArgumentParser(description="Plot the outputs of DNAscent seeBreaks.")
	parser.add_argument("-i", "--input", required=True, help="path to a seeBreaks output file")
	parser.add_argument("-o", "--outdir", default=None, help="directory for the figures (default: alongside the input file)")
	args = parser.parse_args()

	if not os.path.isfile(args.input):
		sys.stderr.write("Error: input file not found: {}\n".format(args.input))
		sys.exit(1)

	outdir = args.outdir if args.outdir else os.path.dirname(os.path.abspath(args.input))
	os.makedirs(outdir, exist_ok=True)
	stem = os.path.splitext(os.path.basename(args.input))[0]

	header, sections = parseSeeBreaks(args.input)

	plotEffect(header, sections, os.path.join(outdir, stem + "_effect.pdf"))
	plotRawFractions(header, sections, os.path.join(outdir, stem + "_rawFractions.pdf"))

	sys.stderr.write("Wrote figures to {}\n".format(outdir))


if __name__ == "__main__":
	main()
