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
names = []
with open('names.txt', 'r') as f:
    for line in f.readlines():
        line = line.replace('_', ' ')  # replace _ with space
        line = line.strip()  # remove newline character
        names.append(line)  # append to list

# DOWNLOAD
sequences = download(names)  # returns dictionary of named sequences

# OUTPUT
counter = 0
for name in sequences.keys():
    seq = sequences[name]  # extract sequence
    if seq:
        seq.name = name  # make sure seq name is name
        seq.id = name  # make sure seq id is name
        seq = seq.format('fasta')  # convert to fasta format
        with open(os.path.join(outdir, "{0}.fasta".format(name)), 'wb') as f:
            f.write("{0}\n".format(seq))
        counter += 1

# FINISH MESSAGE
print("Downloaded [{}] sequences. Finished: [{}]".format(counter, timestamp()))
