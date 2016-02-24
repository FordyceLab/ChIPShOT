import subprocess
import sys
import os
import argparse
from Bio import SeqIO


def index(reference):
    subprocess.call(["bwa", "index", reference])


def fastq_filter(fastq, output, mean_threshold, five_prime_thresh,
                 three_prime_thresh, both_ends_thresh):
    with open(fastq, "r") as fastq_file:
        fq = SeqIO.parse(fastq_file, "fastq")
        for record in fq:
            qualities = record.letter_annotations["phred_quality"]
            mean_quality = sum(qualities)/len(qualities)
            # Have to define mean_threshold
            if mean_quality >= mean_threshold:
                if both_ends_thresh > 0:
                    five_prime_thresh = both_ends_thresh
                    three_prime_thresh = both_ends_thresh
                start = 0
                end = len(record.seq)
                if sum(i >= five_prime_thresh and
                       i >= three_prime_thresh for i in qualities) >= 1:
                    if five_prime_thresh != 0:
                        for i in range(len(qualities)):
                            # Have to define min_base_qual
                            if qualities[i] >= five_prime_thresh:
                                start = i
                                break
                    if three_prime_thresh != 0:
                        for i in range(len(qualities)):
                            if qualities[::-1][i] >= three_prime_thresh:
                                end = - (i + 1)
                                break

                    with open(output, "a") as outfile:
                        SeqIO.write(record[start:end], outfile, "fastq")


def main(argv):
    parser = argparse.ArgumentParser(description="Filter fastq file(s) and \
                                     align to reference genome.")
    parser.add_argument("-f", help="quality cutoff for all bases working from \
                        the 5' end of the read inward", type=int, default=0,
                        nargs="?")
    parser.add_argument("-t", help="quality cutoff for all bases working from \
                        the 3' end of the read inward", type=int, default=0,
                        nargs="?")
    parser.add_argument("-b", help="quality cutoff for all bases working from \
                        both ends of the read inward", type=int, default=0,
                        nargs="?")
    parser.add_argument("-q", help="mean quality cutoff all reads", type=int,
                        default=0, nargs="?")
    parser.add_argument("-r", help="file path to the reference genome",
                        type=str, nargs="?")
    parser.add_argument("-i", help="file path(s) to the fastq files to be \
                        aligned", nargs="*", type=str)
    parser.add_argument("-o", help="path to output file", nargs="?", type=str)

    args = parser.parse_args(argv[1:])

    for fastq in args.i:
        if args.o:
            output = args.o
        else:
            output = fastq.split("\.")[0] + "_filtered.fastq"
        fastq_filter(fastq, output, args.q, args.f, args.t, args.b)

if __name__ == "__main__":
    main(sys.argv)
