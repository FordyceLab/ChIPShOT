import sys
import re
import argparse

# Authored by: Tyler Shimko

def main(argv):
    """
    Clean up the SGD version of the reference genome so that chromosome numbers
    are used.

    Args:
        argv (list) - command line arguments
    """

    # Construct the parser
    parser = argparse.ArgumentParser(description="Download .sra file and \
                                     extract fastq or sam file from GEO")
    parser.add_argument("genome", help="path to reference genome file")
    parser.add_argument("output", help="output file prefix")

    args = parser.parse_args(argv[1:])

    # Open the reference genome and the output
    with open(args.genome, "r") as genome, \
            open(args.output + ".fsa", "w") as output:

        # For each line
        for line in genome:

            # If the line is a separator
            if line.startswith(">"):

                # Get the chromosome information
                chrom = re.search("chromosome=(.+?)]", line)
                if chrom:
                    chrom = re.sub("]", "",
                                   re.sub("chromosome=", "", chrom.group(0)))
                else:
                    print("Chromosome name not found, assuming mito DNA")
                    chrom = "MT"
                line = ">" + chrom + "\n"

            # Write the output file
            output.write(line)

if __name__ == "__main__":
    main(sys.argv)
