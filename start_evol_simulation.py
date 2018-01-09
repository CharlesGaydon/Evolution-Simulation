# coding: utf-8
import argparse

import sys
import evol_simulation as evol

from multiprocessing.pool import Pool
from os import listdir, getcwd
from os.path import isfile, join

# parameters parser


parser = argparse.ArgumentParser(prog='script.py',
                                 usage='%(prog)s files [files ...] --nproc {1,2,3,4,5,6,7,8}',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
Start a simulation using .ini files
-----------------------------------

Simulations can be run with one or more configuration files. Files can
be specified manually in the command line or from a file containing a list 
of config files using the command \'@tasks.txt\'
''',
                                 fromfile_prefix_chars='@')

parser.add_argument('input_files', metavar='files', nargs='+', help='list of config files')
parser.add_argument('--nproc', help='number of processors to use', nargs='?', 
                    type=int, default=4, choices=[1,2,3,4,5,6,7,8], metavar='processors')


#parser.add_argument('N_PROC', metavar='N_PROC', type=int, nargs='+',
                    #help='number of processors to use')
                    
#parser.add_argument('--in', dest='input_file', action='store_const',
                    #const=sum, default=max,
                    #help='input file for the computation')

args = parser.parse_args()

#-------------------------------------------------------------------------------

effective_nproc = min(len(args.input_files), args.nproc)


def PooledComputationOnProcessors(file_list, nproc) :
    
    print('STARTING SIMULATION: %d tasks will be run on %d processors'%(len(file_list), effective_nproc))
    
    results = []
    
    with Pool(effective_nproc) as p :
    
        p.map(evol.start_evol_simulation, file_list)


def PooledComputation(file_list):
    
    return PooledComputationOnProcessors(file_list, effective_nproc)



PooledComputation(args.input_files)




    
    