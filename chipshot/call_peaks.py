import argparse
import subprocess
import sys


def call_peaks(control_bam, treatment_bam, macs_args):
    args = ["macs2", "--nomodel", "-t", treatment_bam]

    if control_bam != "":
        args.extend(["-c", control_bam])

    args.extend(macs_args)

    subprocess.call(args)


def find_pwm_hits():
    return


def recenter_peaks():
    return


def main(argv):
    parser = argparse.ArgumentParser(description="Call ChIP-seq peaks using \
                                     MACS2 and recenter the peaks on the best \
                                     hit from a position weight matrix")
    parser.add_argument("--control", help="control sequencing run, can be \
                        input, IgG pulldown, etc", type=str, nargs="?",
                        metavar="control")
    parser.add_argument("--reference", help="file path to the reference \
                        genome", type=str, nargs="?", metavar="reference")
    parser.add_argument("--extend", help="distance to extend sequence space \
                        on either side of the best hit against the provided \
                        position weight matrix", type=int, default=100,
                        nargs=1, metavar="distance")
    parser.add_argument("input", help="treatment files (ChIP-ed sample) in \
                        bam format", nargs="*", type=str, metavar="bam")
    parser.add_argument("output", help="output file prefix", nargs=1,
                        type=str, metavar="prefix")


if __name__ == "__main__":
    main(sys.argv)
