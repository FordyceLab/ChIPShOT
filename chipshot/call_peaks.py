import argparse
import subprocess
import sys
from Bio import SeqIO, motifs
from Bio.Alphabet import IUPAC
import re
from collections import Counter
import pysam
from tqdm import *


def call_peaks(control_bam, treatment_bam, macs_args, output):
    """
    Provide an interface to MACS2 to call peaks from ChIP-seq data

    Args:
        control_bam (str) - file path to the control bam file
        treatment_bam (str) - file path to the treatment bam file
        macs_args (list) - command line arguments to pass to MACS2
        output (str) - prefix for output files
    """

    # Construct the call
    args = ["macs2", "callpeak", "--nomodel", "-n", output, "-t"]
    args.extend(treatment_bam)

    if control_bam != "":
        args.extend(["-c", control_bam])

    args.extend(macs_args)

    # Call MACS2
    subprocess.call(args)


def calc_coverage(bam):
    bam = pysam.AlignmentFile(bam, "rb")
    coverage = {}
    for line in tqdm(bam):
        try:
            if line.reference_name not in coverage.keys():
                coverage[line.reference_name] = Counter()
            for pos in range(line.reference_start,
                             line.reference_start + len(line.query_sequence)):
                coverage[line.reference_name][pos] += 1
        except:
            pass
    return coverage


def find_pwm_hits(narrow_peak, reference, pfm, output, control_cov, treat_cov):
    """
    Search each peak for the best match against the specified position
    frequency matrix

    Args:
        narrow_peak (str) - path to the narrowPeak file output by MACS2
        reference (str) - file path to the reference genome
        pfm (str) - file path to the position frequency matrix
        output (str) - prefix for the output file
    """

    # Open the peaks and reference genome files
    with open(narrow_peak, "r") as peaks, open(reference, "r") as ref:
        # Parse the reference genome into a dictionary
        records = SeqIO.parse(ref, "fasta", alphabet=IUPAC.unambiguous_dna)
        ref_seq = {record.id: record for record in records}

        # Open and parse the position frequency matrix
        with open(pfm, "r") as pfm:
            matrix = motifs.parse(pfm, "jaspar")[0]
            pwm = matrix.counts.normalize(pseudocounts=.5)
            pssm = pwm.log_odds()

        # Open the output file
        with open(output + "_centeredpeaks.txt", "w") as outfile:

            # Write a line for each centered peak in the output file
            for peak in peaks:
                split_peak = peak.strip().split("\t")
                peak_chrom = split_peak[0]
                peak_start = int(split_peak[1])
                peak_end = int(split_peak[2])
                seq = ref_seq[peak_chrom].seq[peak_start:peak_end]

                hits = [(pos, score) for pos, score in pssm.search(seq)]

                hits.sort(key=lambda hit: hit[1], reverse=True)

                recenter_peak(outfile, ref_seq, peak_chrom, peak_start,
                              peak_end, 100, hits, matrix, control_cov,
                              treat_cov)


def recenter_peak(out_handle, ref_seq, chrom, seq_start, seq_end, slop,
                  hits, matrix, control_cov, treat_cov):
    """
    Recenter peaks on the best hit against the position frequency matrix

    Args:
        out_handle (handle)
        ref_seq (dict) - reference genome sequence stored as a dictionary
        chrom (str) - chromosome identifier
        seq_start (int) - start position of the putative binding sequence
        seq_end (int) - end position of the putative binding sequence
        slop (int) - amount to extend out from the center of the peak on either
            side
        hits (list) - list of hits against pfm in the sequence
        matrix (pfm matrix) - parsed position frequency matrix
    """

    # If hits against the pfm are found
    if hits:

        # Parse the hit for position information
        hit_pos = hits[0][0]
        hit_start = seq_start + hit_pos
        hit_end = seq_start + hit_pos + len(matrix)
        start_adjusted = hit_start - slop
        end_adjusted = hit_end + slop

        # Get the sequence under the peak
        seq = ref_seq[chrom].seq[start_adjusted:end_adjusted]
        rev_seq = seq.reverse_complement()

        # Get the mean coverage across the peak
        coverage = [treat_cov[chrom][pos] - control_cov[chrom][pos] for pos in
                    range(start_adjusted, end_adjusted)]
        mean_cov = sum(coverage) / (end_adjusted - start_adjusted)

        # Determine if hit against psm was on forward or reverse strand
        if hit_pos < 0:
            strand = "-"
            seq = rev_seq
        else:
            strand = "+"

        # Look for the consensus sequence under the peak
        cons_forward = re.search(str(matrix.consensus), str(seq))
        cons_reverse = re.search(str(matrix.consensus), str(rev_seq))
        contains_cons = cons_forward or cons_reverse

        if contains_cons:
            color = "255,0,0"
        else:
            color = "0,0,255"

        # Construct and write the line corresponding to the peak to the output
        # file
        line = [chrom, start_adjusted, end_adjusted, seq, mean_cov, strand, 0,
                0, color]

        line = [str(element) for element in line]
        out_handle.write("\t".join(line) + "\n")


def main(argv):
    """
    Find and recenter peaks around the best hit for a given position weight
    matrix

    Args:
        argv (list) - list of command line arguments
    """

    # Construct an argument parser
    parser = argparse.ArgumentParser(description="Call ChIP-seq peaks using \
                                     MACS2 and recenter the peaks on the best \
                                     hit from a position weight matrix")
    parser.add_argument("--control", help="control sequencing run, can be \
                        input, IgG pulldown, etc", type=str, nargs="?",
                        metavar="control", default="")
    parser.add_argument("--reference", help="file path to the reference \
                        genome", type=str, nargs=1, metavar="reference")
    parser.add_argument("--extend", help="distance to extend sequence space \
                        on either side of the best hit against the provided \
                        position weight matrix", type=int, default=100,
                        nargs=1, metavar="distance")
    parser.add_argument("--pfm", help="position frequency matrix from JASPAR \
                        for the protein of interest", nargs=1, metavar="pfm")
    parser.add_argument("input", help="treatment files (ChIP-ed sample) in \
                        bam format", nargs="*", type=str, metavar="bam")
    parser.add_argument("output", help="output file prefix", nargs=1,
                        type=str, metavar="prefix")

    # Parse command line arguments
    args, macs_args = parser.parse_known_args(argv[1:])

    output = args.output[0]
    reference = args.reference[0]
    pfm = args.pfm[0]

    # Call peaks
    print("Calling ChIP-seq peaks...")
    call_peaks(args.control, args.input, macs_args, output)

    # Get coverage for each position
    print("Calculating coverage for control sample...")
    control_cov = calc_coverage(args.control)
    print("Calculating coverage for control sample...")
    treat_cov = calc_coverage(args.input)

    # Find pfm hits, recenter peaks and write to file
    print("Recentering ChIP-seq peaks...")
    find_pwm_hits(output + "_peaks.narrowPeak", reference, pfm, output,
                  control_cov, treat_cov)

# Read command line arguments if the script is called directly
if __name__ == "__main__":
    main(sys.argv)
