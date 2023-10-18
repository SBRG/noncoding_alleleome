#!/usr/bin/env python3

'''
Created on Sat Aug 06 2022

@author: siddC

Note: This script will not work in its current folder structure.
Put this script in a seperate directory with `prokka_fna_dict.pickle`
present and ensure you have the space needed on your device before
running. Also make sure the location of the ClermonTyping script is
correct.


'''

import os
import pickle
import subprocess

from multiprocessing import Pool


# Location of clermonTyping script
CLERMONT = '/Users/cam/bin/ClermonTyping/clermonTyping.sh'

# Location of genomes
with open('prokka_fna_dict.pickle', 'rb') as f:
    PROKKA_FNA_DICT = pickle.load(f)

assert isinstance(PROKKA_FNA_DICT, dict)


def run_clermonTyping(name, path):
    cmd = f'{CLERMONT} --fasta {path} --name {name}'
    cmd = cmd.split(' ')
    
    subprocess.run(cmd)


if __name__ == '__main__':
    NAME_PATH = [(name, path) for name, path in PROKKA_FNA_DICT.items()]

    with Pool(16) as pool: # number of workers = 16
        pool.starmap(run_clermonTyping, NAME_PATH)

