import sys
import re
import argparse


def main(argv):
    parser = argparse.ArgumentParser(description="Download .sra file and \
                                     extract fastq or sam file from GEO")
    parser.add_argument("genome", help="path to reference genome file")
    parser.add_argument("output", help="output file prefix")

    args = parser.parse_args(argv[1:])

    with open(args.genome, "r") as genome, \
            open(args.output + ".fsa", "w") as output:
        for line in genome:
            if line.startswith(">"):
                chrom = re.search("chromosome=(.+?)]", line)
                if chrom:
                    chrom = re.sub("]", "",
                                   re.sub("chromosome=", "", chrom.group(0)))
                else:
                    print("Chromosome name not found, assuming mito DNA")
                    chrom = "MT"
                line = ">" + chrom

            output.write(line)

if __name__ == "__main__":
    main(sys.argv)
