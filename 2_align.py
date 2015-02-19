#! /bin/usr/env python
# Dom Bennett
# 19/02/2015
# Align sequences

# PACKAGES
import os
from tools import align
from tools import timestamp


# START MESSAGE
print("Align stage, started: [{0}]".format(timestamp()))

# DIRS
indir = '1_download'
outdir = '2_align'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# ALIGN
log = align(infile=os.path.join(indir, 'sequences.fasta'),
            outfile=os.path.join(outdir, 'alignment.fasta'),
            logfile=os.path.join(outdir, 'mafft_log.txt'))
print(log)

# FINISH MESSAGE
print("Finished: [{0}]".format(timestamp()))
