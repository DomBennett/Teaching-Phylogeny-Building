#! /bin/usr/env python
# [YOUR NAME]
# [DATE]
# Tools for phylogeny-building pipeline

from Bio import Entrez  # talk to NCBI's Entrez
from Bio import SeqIO  # read sequences
from random import sample  # random sample function
from datetime import datetime  # for timestamps


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
        term = '"{0}"[PORGN] AND "rbcl"[GENE]'.format(name)
        # search for suitable sequences first
        # create handle to Entrez
        handle = Entrez.esearch(term=term, usehistory='n', retStart=0,
                                retMax=10, retmode="text", db="nucleotide")
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
            handle = Entrez.efetch(id=seqid, rettype='gb', retmode="text",
                                   db="nucleotide")
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
