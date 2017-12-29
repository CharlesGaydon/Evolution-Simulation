# coding: utf-8
import pandas as pd
import simulation
from os import system as oss
import configparser
import time as time
import numpy as np

def import_GFF(gff_file) :
    
    names = ["seqid", "source", "type","start","end","score","strand","phase","attributes"]
    
    seq = pd.read_table(gff_file,sep = '\t', comment='#')
    
    seq_name = seq.columns.values[0]
    seq_length = int(float(seq.columns.values[4]))
    seq.columns = names
    
    return({'seq':seq,'seq_name':seq_name,'seq_length':seq_length})

def load_tab_file(filename):  
    data = pd.read_table(filename, sep='\t', header=0)
    return(data)

def read_config(config_path) :
    
    # Read config 
    
    config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
    # to preserve capital letters
    config.optionxform = str
    config.read(config_path)

    return config

def parse_config(config) :
    
    CONFIG = {}
    
    CONFIG['CONFIG_NAME'] = config.get('INPUTS', 'CONFIG_NAME')
    CONFIG['SIM_TIME'] = config.getint('INPUTS', 'SIM_TIME')
    CONFIG['POP_SIZE'] = config.getint('INPUTS', 'POP_SIZE')
    CONFIG['N_REPS'] = config.getint('INPUTS', 'N_REPS')
    
    CONFIG['INI_PATH'] = config.get('PATHS', 'INI_PATH')
    CONFIG['SIM_PATH'] = config.get('PATHS', 'SIM_PATH')
    CONFIG['ENV_PATH'] = config.get('PATHS', 'ENV_PATH')
    CONFIG['PAR_PATH'] = config.get('PATHS', 'PAR_PATH')
    
    CONFIG['PLASMID_PATH'] = config.get('PATHS', 'PLASMID_PATH')
    CONFIG['TARGET_PATH'] = config.get('PATHS', 'TARGET_PATH')
    
    CONFIG['GFF_FILE'] = config.get('PATHS', 'GFF_FILE')
    CONFIG['TSS_FILE'] = config.get('PATHS', 'TSS_FILE')
    CONFIG['TTS_FILE'] = config.get('PATHS', 'TTS_FILE')
    CONFIG['PROT_FILE'] = config.get('PATHS', 'PROT_FILE')
    
    
    CONFIG['HISTORY_SAVE_NAME'] = config.get('PATHS', 'HISTORY_SAVE_NAME')
    CONFIG['PLASMID_SAVE_NAME'] = config.get('PATHS', 'PLASMID_SAVE_NAME')
    
    CONFIG['PINV'] = config.getfloat('SIMULATION', 'PINV')
    CONFIG['PDEL'] = config.getfloat('SIMULATION', 'PDEL')
    CONFIG['PINS'] = config.getfloat('SIMULATION', 'PINS')
    CONFIG['U'] = config.getint('SIMULATION', 'UNIT')
    CONFIG['ALPHA_C'] = config.getint('SIMULATION', 'ALPHA_C')
    
    CONFIG['STEPS_DONE'] = config.getint('STATE', 'STEPS_DONE')
    CONFIG['REPS_DONE'] = config.getint('STATE', 'REPS_DONE')

    # ---
    
    CONFIG['PROBS'] = np.array([CONFIG['PDEL'],
                                CONFIG['PINS'],
                                CONFIG['PINV']], 
                                dtype=float)
                                
    CONFIG['WPATH'] = CONFIG['SIM_PATH'] + get_working_path(CONFIG) + '/'
    CONFIG['HISTORY_SAVE_PATH'] = CONFIG['WPATH'] + CONFIG['HISTORY_SAVE_NAME']
    CONFIG['PLASMID_SAVE_PATH'] = CONFIG['WPATH'] + CONFIG['PLASMID_SAVE_NAME']

    return CONFIG

# based on simulation.read_config_file_v2
# parse the params.ini file and returns a data object
# for future mnipulation/mutation. 
def read_plasmid_data(config):
    
    # get inputs infos from the config file
    GFF_file = config.get('PATHS', 'GFF_FILE')
    TSS_file = config.get('PATHS', 'TSS_FILE')
    TTS_file = config.get('PATHS', 'TTS_FILE')
    Prot_file = config.get('PATHS', 'PROT_FILE')
    
    TTS = load_tab_file(TTS_file)
    TSS = load_tab_file(TSS_file)
    Prot = load_tab_file(Prot_file)
    GFF = import_GFF(GFF_file)

    data = {'TTS':TTS,'TSS':TSS,'Prot':Prot,'GFF':GFF}
    return(data)

#def make_w_path(probs):
    #p = "pIns_" + str(round(probs[0],2))+"/pDel_" +str(round(probs[1],2))+"/pInv_" +str(round(probs[2],2)) +"/"
    #return(p)

def get_working_path(config) :
    
    probs = config['PROBS']
    alphac = config['ALPHA_C']
    
    #res = '%d_%d_%d_%d'%(round(probs[0],2)*100,
                         #round(probs[1],2)*100,
                         #round(probs[2],2)*100,
                         #alphac)
                         
    res = config['CONFIG_NAME']
    
    return '%s_%s'%(res, time.strftime('%d%m%H%M', time.localtime()))

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
def save_data(data, config):
    
    # Revert old colums names for futher usage of the genome
    names = ['%s'%(data['GFF']['seq']['seqid'][0]), 
             'RefSeq',
             'region',
             '1',
             '%d'%(data['GFF']['seq_length']),
             '.',
             '+',
             '.',
             'ID=id0;Name=%s'%(data['GFF']['seq']['seqid'][0])]
             
    data['GFF']['seq'].columns = names
    
    data["TTS"].to_csv(config['WPATH'] + 'TTS.dat', sep='\t', index=False)
    data["TSS"].to_csv(config['WPATH'] + 'TSS.dat', sep='\t', index=False)
    data["GFF"]["seq"].to_csv(config['WPATH'] + 'gff.gff', sep='\t',index=False)
    data["Prot"].to_csv(config['WPATH'] + 'prot.dat', sep='\t',index=False)
    
    names = ["seqid", "source", "type","start","end","score","strand","phase","attributes"]
    
    data['GFF']['seq'].columns = names