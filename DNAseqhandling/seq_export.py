# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:03:06 2016

@author: Lieneke Kusters
"""

from Bio import SeqIO
import sys
import warnings
import gzip
import shutil
import os


# Functions that can be used to store sequences in a file
def recordwrite(records, filename):
    """ This function will store sequences (in generator format) to a
    possibly zipped (if filename has .gz) fasta (.fa) or fastq (.fastq) file
    and return 1 if succesfull
    """
    cur_version = sys.version_info
    zipit = False

    checkextension = filename.split('.')
    # first check if zipped
    if checkextension[-1] == 'gz':
        zipit = True
        checkextension.pop()  # remove gz
        filename = filename[:-3]

    if checkextension[-1] == 'fna' or checkextension[-1] == 'fa':
        SeqIO.write(records, filename, "fasta")
    elif checkextension[-1] == 'fastq':
        SeqIO.write(records, filename, "fastq")
    else:
        raise ValueError("filename extension {0} not recognised"
                         .format(checkextension[-2]))

    if zipit:
        if cur_version.major == 2:
            # python version is 2, note that I don't expect this to happen!
            warnings.warn("python version is not 3, but {0}"
                          .format(cur_version)
                          )
            with open(filename, 'r') as f_in:
                with gzip.open(filename+'.gz', 'w') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(filename)  # remove the fastq file s.t. only gz remains
        elif cur_version.major == 3:
            with open(filename, 'rb') as f_in:
                with gzip.open(filename+'.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(filename)  # remove the fastq file s.t. only gz remains
    return 1
