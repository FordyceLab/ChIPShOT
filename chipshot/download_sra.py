import subprocess
import sys
import os
import argparse
from urllib.request import urlopen
from Bio import Entrez
import xml.etree.ElementTree as ET


def search(sample_id, email, dump_fq, dump_sam, output):
    """
    Search the Sequence Read Archive for the given sample ID(s) provided by the
    user and download the .sra files for all runs associated with that ID

    Args:
        sample_id (str) - the sample ID to search for in the SRA database
        email (str) - the user's email address (required to search the db)
        dump_fq (bool) - whether to dump the fastq files after downloading the
                         .sra file
        dump_sam (bool) - whether to dump the sam file after downloading the
                          .sra file
    """

    if email == "":
        email = input("Please enter your email to search the GEO database: ")

    # Set the emial
    Entrez.email = email

    # Print the sample ID being searched
    print('Searching for sample ID "{}"'.format(sample_id))

    # Search the Sequence Read Archive database for the sample ID
    handle = Entrez.esearch(db="sra", term=sample_id)

    # Get the ID number for the sample
    ids = Entrez.read(handle)['IdList']

    # Exit if the sample search has no hits or more than one hit
    if len(ids) == 0:
        print('No records matching ID search term: "{}"'.format(sample_id))
        return
    if len(ids) > 1:
        print('ID search term "{}" is ambiguous'.format(sample_id))
        return

    # Retrieve further information about the sample search
    handle = Entrez.efetch(db="sra", id=ids[0])

    # Parse the returned XML
    tree = ET.parse(handle)
    root = tree.getroot()
    runs = root.find("EXPERIMENT_PACKAGE").find("RUN_SET").findall("RUN")

    # Print a warning if the sample has multiple sequencing runs
    if len(runs) > 1:
        print('Multiple runs found for sample ID "{}"'.format(sample_id))

    # Make the output directory
    dirname = "./{}".format(output)

    if not os.path.exists(dirname):
        os.mkdir(dirname)

    # Get the individual run accession number
    accessions = [run.get("accession") for run in runs]

    if len(accessions) > 1:
        append_acc = True
    else:
        append_acc = False

    # For each accession number in the list of runs
    for accession in accessions:
        print('Downloading sequence read archive for accession "{acc}"'
              .format(acc=accession))

        if append_acc:
            f = output + "_" + accession
        else:
            f = output

        # Download the .sra file and save in the approppriate directory
        filename = "./{}/{}.sra".format(output, f)
        download(accession, filename)

        # If the user wanted to dump the fastq files or the sam file, do so
        if dump_fq:
            subprocess.call(["fastq-dump", "--split-files", "-O",  dirname,
                             filename])

        if dump_sam:
            subprocess.call(["sam-dump", "--output-file",
                             "{}/{}.sam".format(dirname, f), filename])


def download(accession, outfile):
    """
    Subroutine to download the specific run as a .sra file

    Args:
        accession (str) - accession number for the run
        outfile (str) - the name of the file to write
    """

    # Construct the url
    url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"
    acc = "{first3}/{first6}/{full}/{full}.sra".format(first3=accession[:3],
                                                       first6=accession[:6],
                                                       full=accession)
    url = url + acc

    # Download and save the .sra file
    with open(outfile, "wb") as output:
        response = urlopen(url)
        output.write(response.read())


def main(argv):
    """
    Run the search and download subroutines for the sample ID specified by the
    user in the command line areguments

    Args:
        argv (list) - list of command line arguments
    """

    # Construct a parser and parse all of the command line arguments
    parser = argparse.ArgumentParser(description="Download .sra file and \
                                     extract fastq or sam file from GEO")
    parser.add_argument("--fastq", help="dump fastq files",
                        action="store_true")
    parser.add_argument("--sam", help="dump sam file", action="store_true")
    parser.add_argument("--email", help="email address; required by GEO to \
                        search their database", type=str, nargs=1,
                        metavar="email")
    parser.add_argument("ids", help="sample ID numbers from GEO", nargs="*",
                        type=str, metavar="ID")
    parser.add_argument("output", help="output file prefix", nargs=1,
                        metavar="prefix")

    args = parser.parse_args(argv[1:])
    email = args.email[0]

    # Print an error and exit if a sample ID is not provided
    if not args.ids:
        print("Please provide a sample ID")
        sys.exit(1)

    # For each sample ID, search and download
    for sample_id in args.ids:
        search(sample_id, email, args.fastq, args.sam, args.output[0])

# Read command line arguments if the script is called directly
if __name__ == "__main__":
    main(sys.argv)
