import subprocess
import sys
import os
import argparse
from Bio import SeqIO


def fastq_filter(fastq, output, mean_threshold, five_prime_thresh,
                 three_prime_thresh, both_ends_thresh):
    if (mean_threshold != 0 or five_prime_thresh != 0 or
            three_prime_thresh != 0 or both_ends_thresh != 0):
        with open(fastq, "r") as fastq_file, open(output, "a") as outfile:
            fq = SeqIO.parse(fastq_file, "fastq")
            for record in fq:
                qualities = record.letter_annotations["phred_quality"]
                mean_quality = sum(qualities)/len(qualities)
                if mean_quality >= mean_threshold:
                    if both_ends_thresh > 0:
                        five_prime_thresh = both_ends_thresh
                        three_prime_thresh = both_ends_thresh
                    start = 0
                    end = len(record.seq)
                    if sum(i >= five_prime_thresh and
                           i >= three_prime_thresh for i in qualities) > 0:
                        if five_prime_thresh != 0:
                            for i in range(len(qualities)):
                                if qualities[i] >= five_prime_thresh:
                                    start = i
                                    break
                        if three_prime_thresh != 0:
                            for i in range(len(qualities)):
                                if qualities[::-1][i] >= three_prime_thresh:
                                    end = - (i + 1)
                                    break

                    SeqIO.write(record[start:end], outfile, "fastq")
        return output
    else:
        return fastq


def index(reference):
    subprocess.call(["bwa", "index", reference])


def align(reference, fastq, output, bwa_args):
    index(reference)
    args = ["bwa", "mem"]
    args.extend(bwa_args)
    args.append(reference)
    args.extend(fastq)
    print(args)

    with open(output, "w") as outfile:
        subprocess.call(args, stdout=outfile)


def main(argv):
    parser = argparse.ArgumentParser(description="Filter fastq file(s) and \
                                     align to reference genome.")
    parser.add_argument("--five", help="quality cutoff for all bases working from \
                        the 5' end of the read inward", type=int, default=0,
                        nargs="?", metavar="quality")
    parser.add_argument("--three", help="quality cutoff for all bases working from \
                        the 3' end of the read inward", type=int, default=0,
                        nargs="?", metavar="quality")
    parser.add_argument("--both", help="quality cutoff for all bases working from \
                        both ends of the read inward", type=int, default=0,
                        nargs="?", metavar="quality")
    parser.add_argument("--quality", help="mean quality cutoff all reads",
                        type=int, default=0, nargs="?", metavar="quality")
    parser.add_argument("--min-length", help="minimum length of reads after \
                        all trimming is performed", type=int, default=0,
                        nargs="?", metavar="length")
    parser.add_argument("--reference", help="file path to the reference \
                        genome", type=str, nargs="?", metavar="reference")
    parser.add_argument("input", help="file path(s) to the fastq files to be \
                        aligned", nargs="*", type=str, metavar="fastq")
    parser.add_argument("output", help="output file prefix", nargs=1,
                        type=str, metavar="prefix")

    args, bwa_args = parser.parse_known_args(argv[1:])

    output_prefix = args.output[0]
    output = output_prefix + ".sam"

    pair = 1

    filtered_fastqs = []

    for fastq in args.input:
        if len(args.input) == 1:
            fastq_output = output_prefix + ".fastq"
        else:
            fastq_output = "{}-{}.fastq".format(output_prefix, pair)

        fastq_file = fastq_filter(fastq, fastq_output, args.quality, args.five,
                                  args.three, args.both)

        filtered_fastqs.append(fastq_file)

        pair += 1

    if args.reference:
        align(args.reference, filtered_fastqs, output, bwa_args)

if __name__ == "__main__":
    main(sys.argv)
