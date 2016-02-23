import subprocess
import sys
import os
import getopt
from urllib.request import urlopen
from Bio import Entrez
import xml.etree.ElementTree as ET


def search(gsm, email, dump_fq, dump_sam):
    """
    Search the Sequence Read Archive for the given GSM ID(s) provided by the
    user and download the .sra files for all runs associated with that ID

    Args:
        gsm (str) - the GSM id to search for in the SRA database
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

    # Print the GSM ID being searched
    print('Searching for GSM "{gsm}"'.format(gsm=gsm))

    # Search the Sequence Read Archive database for the GSM ID
    handle = Entrez.esearch(db="sra", term=gsm)

    # Get the ID number for the GSM
    ids = Entrez.read(handle)['IdList']

    # Exit if the GSM search has no hits or more than one hit
    if len(ids) == 0:
        print('No records matching GSM search term: "{gsm}"'.format(gsm=gsm))
        return
    if len(ids) > 1:
        print('GSM search term "{gsm}" is ambiguous'.format(gsm=gsm))
        return

    # Retrieve further information about the GSM search
    handle = Entrez.efetch(db="sra", id=ids[0])

    # Parse the returned XML
    tree = ET.parse(handle)
    root = tree.getroot()
    runs = root.find("EXPERIMENT_PACKAGE").find("RUN_SET").findall("RUN")

    # Print a warning if the GSM has multiple sequencing runs
    if len(runs) > 1:
        print('Multiple runs found for GSM ID "{gsm}"'.format(gsm=gsm))

    # Make the output directory
    dirname = "./{gsm}".format(gsm=gsm)

    if not os.path.exists(dirname):
        os.mkdir(dirname)

    # Get the individual run accession number
    accessions = [run.get("accession") for run in runs]

    # For each accession number in the list of runs
    for accession in accessions:
        print('Downloading sequence read archive for accession "{acc}"'
              .format(acc=accession))

        # Download the .sra file and save in the approppriate directory
        filename = "./{gsm}/{acc}.sra".format(gsm=gsm, acc=accession)
        download(accession, filename)

        # If the user wanted to dump the fastq files or the sam file, do so
        if dump_fq:
            subprocess.call(["fastq-dump", "--split-files", "-O",  dirname,
                             filename])

        if dump_sam:
            subprocess.call(["sam-dump", "--output-file",
                             "{}/{}.sam".format(dirname, accession), filename])


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
    Run the search and download subroutines for the GSM ID specified by the
    user in the command line areguments

    Args:
        argv (list) - list of command line arguments
    """

    # Set the default arguments
    dump_fq = False
    dump_sam = False
    email = ""

    # Get command line arguments if none are supplied
    if argv == "":
        argv = sys.argv
    else:
        argv = argv[1:]

    # Parse the arguments
    opts, args = getopt.getopt(argv, shortopts="e:fs")

    # Set the arguments based on the command line parameters
    for opt, arg in opts:
        if opt == "-f":
            dump_fq = True
        if opt == "-s":
            dump_sam = True
        if opt == "-e":
            email = arg

    # For each GSM ID, search and download
    for arg in args:
        search(arg, email, dump_fq, dump_sam)


if __name__ == "__main__":
    main(sys.argv)
