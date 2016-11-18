# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:03:06 2016

@author: Lieneke Kusters
"""

import gzip
from Bio import SeqIO
import sys
import warnings


# Functions that can be used to load sequences from a file
def seqgenerator(filenames_list):
    """ This function will load sequences from a zipped fasta or fastq file
    and return a generator of sequences (without identifiers)
    """
    cur_version = sys.version_info

    if cur_version.major == 2:
        # python version is 2, note that I don't expect this to happen!
        warnings.warn("python version is not 3, but {0}"
                      .format(cur_version)
                      )
        for filename in filenames_list:
            handle, fileformat = _openile2(filename)
            for record in SeqIO.parse(handle, fileformat):
                # sys.stderr.write(str(record.description)+'\n')
                yield str(record.seq)

    elif cur_version.major == 3:
        for filename in filenames_list:
            handle, fileformat = _openile(filename)
            for record in SeqIO.parse(handle, fileformat):
                # sys.stderr.write(str(record.description)+'\n')
                yield str(record.seq)

    else:
        raise EnvironmentError("python version {0} incompatible with code"
                               .format(cur_version))


def recordgenerator(filenames_list):
    """ This function will load sequences from a zipped fasta or fastq file
    and return a generator of sequence records
    """
    cur_version = sys.version_info

    if cur_version.major == 2 and cur_version.minor == 7:
        for filename in filenames_list:
            handle, fileformat = _openile2(filename)
            for record in SeqIO.parse(handle, 'fasta'):
                yield record

    elif cur_version.major == 3:
        for filename in filenames_list:
            handle, fileformat = _openile(filename)
            for record in SeqIO.parse(handle, 'fasta'):
                yield record

    else:
        raise EnvironmentError("python version {0} incompatible with code"
                               .format(cur_version))


def _openile(filename):
    """ Open file (that is possibly zipped) and return file handle and
    extension

    We assume that the version is python3
    1) we check whether it is zip-file and unzip if needbe
    2) we check fasta or fastq file format
    """

    checkextension = filename.split('.')
    if checkextension[-1] == 'gz':
        handle = gzip.open(filename, 'rt')
        checkextension.pop(-1)
    else:
        handle = open(filename, 'rt')

    if checkextension[-1] == 'fna' or checkextension[-1] == 'fa':
        fileformat = 'fasta'
    elif checkextension[-1] == 'fastq':
        fileformat = 'fastq'
    else:
        raise ValueError("filename extension {0} not recognised"
                         .format(checkextension[-1]))
    return handle, fileformat


def _openile2(filename):
    """ Open file for python2
    Open file (that is possibly zipped) and return file handle and
    extension

    We assume that the version is python3
    1) we check whether it is zip-file and unzip if needbe
    2) we check fasta or fastq file format
    """

    checkextension = filename.split('.')
    if checkextension[-1] == 'gz':
        handle = gzip.open(filename, 'r')
        checkextension.pop(-1)
    else:
        handle = open(filename, 'r')

    if checkextension[-1] == 'fna' or checkextension[-1] == 'fa':
        fileformat = 'fasta'
    elif checkextension[-1] == 'fastq':
        fileformat = 'fastq'
    else:
        raise ValueError("filename extension {0} not recognised"
                         .format(checkextension[-1]))
    return handle, fileformat
