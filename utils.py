# coding: utf-8
import pandas as pd
import simulation
from os import system as oss
import configparser


def import_GFF(gff_file) :
    seq = pd.read_table(gff_file,sep = '\t', comment='#')
    seq_name = seq.columns.values[0]
    seq_length = int(seq.columns.values[4])
    return({'seq':seq,'seq_name':seq_name,'seq_length':seq_length})

def load_tab_file(filename):  
    data = pd.read_table(filename, sep='\t', header=0)
    return(data)

# based on simulation.read_config_file_v2
# parse the params.ini file and returns a data object
# for future mnipulation/mutation. 
def import_data_from_params_seq_file(path):

    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    config.read(path)
    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')
    TTS = load_tab_file(TTS_file)
    TSS = load_tab_file(TSS_file)
    Prot = load_tab_file(Prot_file)
    GFF = import_GFF(GFF_file)

    data = {'TTS':TTS,'TSS':TSS,'Prot':Prot,'GFF':GFF}
    return(data)

def copy_to_working_path(path_params_seq, working_path):
    #modify the config file from path_params_seq and copy it in working_path.
    config = simulation.read_config_file(path_params_seq+"params_seq.ini")
    config.set('INPUTS','TSS', working_path+"TSS.dat")
    config.set('INPUTS','TTS', working_path+"TTS.dat")
    config.set('INPUTS','GFF', working_path+"gff.gff")
    config.set('INPUTS','BARR_FIX', working_path+"prot.dat")
    #write it in working_path
    with open(working_path+"params_seq.ini", 'w') as configfile:
        config.write(configfile)
    #copy the files.
    oss("cp "+path_params_seq+"TSS.dat "+working_path+"TSS.dat")
    oss("cp "+path_params_seq+"TTS.dat "+working_path+"TTS.dat")
    oss("cp "+path_params_seq+"prot.dat "+working_path+"prot.dat")
    oss("cp "+path_params_seq+"gff.gff "+working_path+"gff.gff")

# Removed 'h=True' because this argument is not known by my Pandas version ...
# Seems to work
def save_data_to_path(data,working_path):
    data["TTS"].to_csv(working_path+"TTS.dat", sep='\t', index=False)
    data["TSS"].to_csv(working_path+"TSS.dat", sep='\t', index=False)
    data["GFF"]["seq"].to_csv(working_path+"gff.gff", sep='\t',index=False)
    data["Prot"].to_csv(working_path+"prot.dat", sep='\t',index=False)