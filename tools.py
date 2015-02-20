#! /bin/usr/env python
# [Your name]
# [Date]
'''
Tools for phylogeny-building pipeline
'''

# PACKAGES
import os
import re
from Bio import Entrez  # talk to NCBI's Entrez
from Bio import SeqIO  # read sequences
from random import sample  # random sample function
from random import randint  # random number generator
from datetime import datetime  # for timestamps
from Bio.Align.Applications import MafftCommandline  # MAFFT command-line
from Bio.Phylo.Applications import RaxmlCommandline  # RAxML command-line


# FUNCTIONS
def timestamp():
    '''Return formatted timestamp'''
    timestamp = datetime.today().strftime("%A, %d %B %Y %I:%M%p")
    return timestamp


def download(names):
    '''Return dictionary of names and sequences'''
    sequences = {}
    for name in names:
        # print progress
        print('.... downloading for [{0}]'.format(name))
        # build a search term: consisting of principle organism and gene
        # for more fields see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
        # TODO -- create search term for name and rbcl
        # search for suitable sequences first
        # create handle to Entrez
        # TODO -- run an Entrez.esearch function with term
        # read handle
        results = Entrez.read(handle)
        # always close handle
        handle.close()
        # extract suitable seqids from results dict
        seqids = results['IdList']
        # choose one at random
        if seqids:
            seqid = sample(seqids, 1)[0]
            # fetch the sequence
            # TODO -- run an Entrez.efetch to pull down sequences
            # rettype is gb, must parse like this
            results_iter = SeqIO.parse(handle, 'gb')
            result = [x for x in results_iter][0]
            handle.close()
        else:
            print('........ no sequences found')
            result = None
        # add to dictionary of name and sequence
        sequences[name] = result
    return sequences


def align(infile, outfile, logfile):
    '''Write MAFFT generated alignment to outfile from sequences in infile.
    Return MAFFT log.'''
    # use 'which' to find your mafft
    # generate a command-line
    # TODO -- run MAFFT commandline function
    # print the commandline that will be used
    print 'Using command-line:\n.... {0}'.format(cline)
    # run command-line
    stdout, stderr = cline()
    # write-out results
    with open(outfile, 'w') as f:
        f.write(stdout)
    with open(logfile, 'w') as f:
        f.write(stderr)
    return stderr


def estimate(infile, outfile, logfile):
    '''Write RAxML generated tree to outfile from alignment in infile.
    Return RAxML log.'''
    # TODO -- run RAxML commandline function
    # print the commandline that will be used
    print 'Using command-line:\n.... {0}'.format(cline)
    # run command-line
    stdout, stderr = cline()
    # move best tree to outfile
    os.rename('RAxML_bestTree.temp', outfile)
    # remove raxml outfiles
    files = os.listdir(os.getcwd())
    for file in files:
        if re.search("(RAxML)", file):
            os.remove(file)
        if re.search("\.reduced$", file):
            os.remove(file)
    # write-out
    # TODO -- write out and return log
