#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import os


def gen_bins(fastafile, resultfile, outputdir):
    # read fasta file
    print("Processing file:\t{}".format(fastafile))
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    print("Reading Map:\t{}".format(resultfile))
    dic = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name] = []
                dic[cluster_name].append(contig_name)
    print("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    # Number bins sequentially (0, 1, 2, ...), one number per cluster/bin — NOT
    # per contig (the old code incremented inside the contig loop, producing
    # cumulative-contig-count IDs and empty files for all-missing clusters).
    bin_name = 0
    for _, cluster in dic.items():
        lines = []
        for contig_name in cluster:
            sequence = sequences.get(">" + contig_name)
            if sequence is None:
                continue
            lines.append(">" + contig_name + "\n")
            lines.append(sequence + "\n")
        if not lines:
            continue
        binfile = os.path.join(outputdir, "{}.fa".format(bin_name))
        with open(binfile, "w") as f:
            f.writelines(lines)
        bin_name += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="original fasta file")
    parser.add_argument("-r", help="tsv version result file")
    parser.add_argument("-o", help="output dir")
    args = parser.parse_args()
    gen_bins(args.f, args.r, args.o)
