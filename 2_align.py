#! /bin/usr/env python
# [Your name]
# [Date]
'''
Aligns sequences
'''

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
# TODO -- run align()

# FINISH MESSAGE
print("Finished: [{0}]".format(timestamp()))
