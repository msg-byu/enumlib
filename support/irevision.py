import re
import argparse

rx = r'(?P<prefix><fortpy[\s\w"=\d]+revision=)"(\d+)"'

def get_revision():
    with open(args["revision"]) as f:
        revision = re.match("(?P<revision>\d+)\w", f.read().strip()).group("revision")

    return int(revision)

def initialize():
    revision = get_revision()

    with open(args["file"]) as f:
        contents = f.read()
        contents = re.sub(rx, '\g<prefix>"{}"'.format(revision), contents)

    with open(args["file"], 'w') as f:
        f.write(contents)

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy File Comparison Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("file", help="Specify the path to the file to update revision numbers for.")
parser.add_argument("revision", help="Specify the file name that contains the revision number.")
#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())
initialize()
