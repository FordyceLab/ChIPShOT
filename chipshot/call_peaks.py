import argparse
import subprocess
import sys
from Bio import SeqIO, motifs
from Bio.Alphabet import IUPAC
import re


def call_peaks(control_bam, treatment_bam, macs_args, output):
    args = ["macs2", "callpeak", "--nomodel", "-n", output, "-t"]
    args.extend(treatment_bam)

    if control_bam != "":
        args.extend(["-c", control_bam])

    args.extend(macs_args)

    subprocess.call(args)


def find_pwm_hits(narrow_peak, reference, pfm, output):
    with open(narrow_peak, "r") as peaks, open(reference, "r") as ref:
        records = SeqIO.parse(ref, "fasta", alphabet=IUPAC.unambiguous_dna)
        ref_seq = {record.id: record for record in records}
        with open(pfm, "r") as pfm:
            matrix = motifs.parse(pfm, "jaspar")[0]
            pwm = matrix.counts.normalize(pseudocounts=.5)
            pssm = pwm.log_odds()

        with open(output + "_centeredpeaks.txt", "w") as outfile:
            header = ["CHROM", "START", "END", "STRAND", "CONTAINS_CONSENSUS",
                      "HIT_SEQ", "SEQ"]
            outfile.write("\t".join(header) + "\n")
            for peak in peaks:
                split_peak = peak.strip().split("\t")
                peak_chrom = split_peak[0]
                peak_start = int(split_peak[1])
                peak_end = int(split_peak[2])
                seq = ref_seq[peak_chrom].seq[peak_start:peak_end]

                hits = [(pos, score) for pos, score in pssm.search(seq)]

                hits.sort(key=lambda hit: hit[1], reverse=True)

                recenter_peak(outfile, ref_seq, peak_chrom, peak_start,
                              peak_end, 100, hits, matrix)


def recenter_peak(out_handle, ref_seq, chrom, seq_start, seq_end, slop,
                  hits, matrix):
    if hits:
        hit_pos = hits[0][0]
        hit_start = seq_start + hit_pos
        hit_end = seq_start + hit_pos + len(matrix)
        start_adjusted = hit_start - slop
        end_adjusted = hit_end + slop

        seq = ref_seq[chrom].seq[start_adjusted:end_adjusted]
        rev_seq = seq.reverse_complement()
        cons_forward = re.search(str(matrix.consensus), str(seq))
        cons_reverse = re.search(str(matrix.consensus), str(rev_seq))
        if cons_forward or cons_reverse:
            contains_cons = 1
        else:
            contains_cons = 0

        if hit_pos < 0:
            strand = "-"
            seq = rev_seq
        else:
            strand = "+"

        hit_seq = seq[slop:slop + len(matrix)]

        line = [chrom, start_adjusted, end_adjusted, strand, contains_cons,
                hit_seq, seq]

        line = [str(element) for element in line]
        out_handle.write("\t".join(line) + "\n")


def main(argv):
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

    args, macs_args = parser.parse_known_args(argv[1:])

    output = args.output[0]
    reference = args.reference[0]
    pfm = args.pfm[0]

    call_peaks(args.control, args.input, macs_args, output)

    find_pwm_hits(output + "_peaks.narrowPeak", reference, pfm, output)


if __name__ == "__main__":
    main(sys.argv)
