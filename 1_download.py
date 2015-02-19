#! /bin/usr/env python
# [YOUR NAME]
# [DATE]
# Download sequences from GenBank

# PACKAGES
import os  # create, control and modify dirs and files
from Bio import Entrez  # talk to NCBI's Entrez
from tools import download  # read in your own functions like this
from tools import timestamp

# START MESSAGE
print("Download stage, started: [{}]".format(timestamp()))

# DIRS
outdir = '1_download'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# PARAMETERS
Entrez.email = 'your.email@address.here'

# INPUT
with open('names.txt', 'r') as f:
    # read each line, set as new item in list, remove \n
    names = f.read().splitlines()

# DOWNLOAD
sequences = download(names)  # returns dictionary of named sequences

# OUTPUT
counter = 0
with open(os.path.join(outdir, "sequences.fasta"), 'wb') as f:
    for name in sequences.keys():
        seq = sequences[name]  # extract sequence
        if seq:
            seq.name = name  # make sure seq name is name
            seq.id = name  # make sure seq id is name
            seq = seq.format('fasta')  # convert to fasta format
            f.write("{0}\n".format(seq))
            counter += 1

# FINISH MESSAGE
print("Downloaded [{}] sequences. Finished: [{}]".format(counter, timestamp()))
