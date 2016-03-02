import subprocess
import sys
import argparse
from Bio import SeqIO
from tqdm import *


def fastq_filter(fastq, output, mean_threshold, five_prime_thresh,
                 three_prime_thresh, both_ends_thresh, min_length):
    """
    Filter a given fastq file based on average quality and minimum length, as
    well as provide functionality to trim reads from 5' or 3' ends until a base
    of minimum quality is reached

    Args:
        fastq (str) - path to the fastq file of interest
        output (str) - path to output (filtered) fastq file
        mean_threshold (int) - minimum average quality
        five_prime_thresh (int) - minimum base quality to trim to from the 5'
            end of the read
        three_prime_thresh (int) - minimum base quality to trim to from the 3'
            end of the read
        both_ends_thresh (int) - minimum base quality to trim to from both
            ends of the read
        min_length (int) - minimum lenth of the read after all trimming

    Returns:
        The name of the output file
    """

    # Ensure that the some filtering parameter is set
    if (mean_threshold != 0 or five_prime_thresh != 0 or
            three_prime_thresh != 0 or both_ends_thresh != 0):

        # Open the input and fastq file as well the output (filtered) file
        with open(fastq, "r") as fastq_file, open(output, "a") as outfile:

            # Parse the fastq file
            fq = SeqIO.parse(fastq_file, "fastq")

            # Iterate through the records with a nice little progress meter
            for record in tqdm(fq):

                # Get the base qualities and calculate the mean quality score
                # for the read
                qualities = record.letter_annotations["phred_quality"]
                mean_quality = sum(qualities)/len(qualities)

                # Check if the mean quality exceeds the threshold
                if mean_quality >= mean_threshold:

                    # If both ends threshold is set, it trumps either
                    # individual end threshold
                    if both_ends_thresh > 0:
                        five_prime_thresh = both_ends_thresh
                        three_prime_thresh = both_ends_thresh

                    # Define the original start and end of the read before
                    # trimming
                    start = 0
                    end = len(record.seq)

                    # Make sure the read has a base with a quality score over
                    # both of the trimming thresholds
                    if sum(i >= five_prime_thresh and
                           i >= three_prime_thresh for i in qualities) > 0:

                        # Trim back from the 5' end
                        if five_prime_thresh != 0:
                            for i in range(len(qualities)):
                                if qualities[i] >= five_prime_thresh:
                                    start = i
                                    break

                        # Trim back from the 3' end
                        if three_prime_thresh != 0:
                            for i in range(len(qualities)):
                                if qualities[::-1][i] >= three_prime_thresh:
                                    end = - (i + 1)
                                    break

                    # If the trimmed read is longer than the minimum length,
                    # write that sequence to the output file
                    if len(record[start:end]) >= min_length:
                        SeqIO.write(record[start:end], outfile, "fastq")

        # Return the name of the output file (if any filtering was performed)
        return output
    else:
        # Return the name of the original fastq (if no filtering was performed)
        return fastq


def index(reference):
    """
    Call `bwa index` on the selected reference genome sequence

    Args:
        reference (str) - path to reference genome fasta file
    """

    subprocess.call(["bwa", "index", reference])


def align(reference, fastq, output, bwa_args):
    """
    Align a fastq file to a given reference genome

    Args:
        reference (str) - path to the reference genome
        fastq (str) - path to the fastq file to align
        output (str) - prefix for the output bam file
        bwa_args (list) - list of arguments to pass thorough to bwa
    """

    # Index the reference genome
    index(reference)

    # Construct the argument list to pass to bwa
    args = ["bwa", "mem"]
    args.extend(bwa_args)
    args.append(reference)
    args.extend(fastq)

    # Call bwa and pipe output to samtools for sam to bam conversion and write
    # to output file
    with open("{}.bam".format(output), "wb") as outfile:
        alignment = subprocess.Popen(args, stdout=subprocess.PIPE)
        subprocess.call(["samtools", "view", "-bS", "-"],
                        stdin=alignment.stdout, stdout=outfile)


def main(argv):
    """
    Run the search and download subroutines for the GSM ID specified by the
    user in the command line areguments

    Args:
        argv (list) - list of command line arguments
    """

    # Construct the parser 
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
    parser.add_argument("--minlength", help="minimum length of reads after \
                        all trimming is performed", type=int, default=0,
                        nargs="?", metavar="length")
    parser.add_argument("--reference", help="file path to the reference \
                        genome", type=str, nargs="?", metavar="reference")
    parser.add_argument("input", help="file path(s) to the fastq files to be \
                        aligned", nargs="*", type=str, metavar="fastq")
    parser.add_argument("output", help="output file prefix", nargs=1,
                        type=str, metavar="prefix")

    # Parse the command line arguments
    args, bwa_args = parser.parse_known_args(argv[1:])

    output_prefix = args.output[0]

    # Set up to handle paired reads in two files
    pair = 1

    filtered_fastqs = []

    # For each fastq given (can be muliptle if paired end data)
    for fastq in args.input:
        print("Filtering " + fastq)

        # Handle naming of output file for one or two input files
        if len(args.input) == 1:
            fastq_output = output_prefix + ".fastq"
        else:
            fastq_output = "{}-{}.fastq".format(output_prefix, pair)

        # Filter the fastq
        fastq_file = fastq_filter(fastq, fastq_output, args.quality, args.five,
                                  args.three, args.both, args.minlength)

        # Make a list of the filtered fastq file paths
        filtered_fastqs.append(fastq_file)

        pair += 1

    # If a reference genome is included as an argument, perform an alignment
    if args.reference:
        align(args.reference, filtered_fastqs, output_prefix, bwa_args)

# Read command line arguments if the script is called directly
if __name__ == "__main__":
    main(sys.argv)
