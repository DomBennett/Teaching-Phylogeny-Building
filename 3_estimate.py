#! /bin/usr/env python
# Dom Bennett
# 19/02/2015
'''
Generates trees with RAxML
'''

# PACKAGES
import os
from tools import estimate
from tools import timestamp

# START MESSAGE
print("Estimate stage, started: [{0}]".format(timestamp()))

# DIRS
indir = '2_align'
outdir = '3_estimate'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# ESTIMATE
log = estimate(infile=os.path.join(indir, 'alignment.fasta'),
               outfile=os.path.join(outdir, 'best_tree.tre'),
               logfile=os.path.join(outdir, 'RAxML_log.txt'))
print(log)

# FINISH MESSAGE
print("Finished: [{0}]".format(timestamp()))
