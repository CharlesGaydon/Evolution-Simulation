# coding: utf-8
import pandas as pd
import numpy as np
import copy
import os
import simulation
import utils
PARAMS = {}

def start_evol_simulation(INI_file) :

    #import model parameters from INI_fileparams_evol.init (dont le nom du fichier params.init)
    #note : utiliser syntaxe de simulation.read_config_file(path) ?
    PARAMS['U'] = 150
    PARAMS['POP_SIZE'] = 1 #pour l'instant on ne gere que ce cas. N>1 plus tard.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  [0.33,0.33,0.33]
    PARAMS['w_path_1'] = "tousgenesidentiques/N_1/pIns_0.33/pDel_0.33/pInv_0.33/1/" #TODO : automatiser
    PARAMS['w_path_2'] = "tousgenesidentiques/N_1/pIns_0.33/pDel_0.33/pInv_0.33/2/"
    os.system("mkdir -p " + PARAMS['w_path_1']) #à adapter et peut etre deplacer dans plasmid.
    os.system("mkdir -p " + PARAMS['w_path_2'])

    Plasmid = plasmid()
    #créer une population et la stocker
        #chaque individu est crée avec le fichier de params.ini comme argument, #PAS ENCORE car N=1

    #boucler SIM_TIME fois (en nombre de mut ?)
        #tirer les mutations possibles et les exécuter ; 
    Plasmid.mutate()
        #lancer la simulation : mutation, (N=1 : création d'un fichier next_params.ini et des fichiers 
        #TTS TSS GFF Prot ... _next correspondant)
        #N!=1 : actualisation des fichiers de donné.
        # le calcul de la fitness sur params_next, sa comparaison,  au sein d'un individu. 
        #Note : faire le tri dans les outputs enregistrés de simulation.py
        #actualiser ensuite les individus (si un seul individu inclure Metropolis dans l'individu.)
        #si plusieurs individus, reprod sexuelle impossible... duplication simple ?
        #differencier les logs si N = 1 ou N > 1 
    #ptit plot de l'évolution de la fitness max qui fait plaiz'



class plasmid:

    def __init__(self, ID = 1):
        self.ID = ID
        self.probs = PARAMS['probs']
        self.w_path_1 = PARAMS['w_path_1']
        self.w_path_2 = PARAMS['w_path_2']
        self.o_path = PARAMS['path_params_seq']
        # self.path_params_seq = PARAMS['path_params_seq'] #useless ?
        # copy of files in the working paths
        utils.copy_to_working_path(PARAMS['path_params_seq'], self.w_path_1)
        utils.copy_to_working_path(PARAMS['path_params_seq'], self.w_path_2)
        self.data = utils.import_data_from_params_seq_file(self.o_path+"params_seq.ini")
        
        
        #lancer une première simulation pour calculer fitness de w_path_1
        self.fitness = 0
        self.hist_fitness = [self.fitness]
    def save(self) :
        #to save hist_fitness, and maybe other stuff, later
        pass

    def mutate(self) : 
        #tirage et choix de la mut
        #  generation des updated data et
        # sauvegarde dans files de w_path_2
        # lancer simulation avec params_seq de w_path_2
        
        #test
        print(self.data)
        updated_data = self.U_deletion()
        print(updated_data)
        utils.save_data_to_path(updated_data,self.w_path_2)

        #calculer fitness
        #comparer et conserver w1 ou w2
        #actualiser en conséquence.

        # puis
        # enregistrer la fitness dans self.fitness et updater l'historique
        # enfin
        # update le self.data choisi for next mutation
        pass

    def U_deletion(self) :
        data = self.data
        l = data['GFF']['seq_length']
        #localisation
        start=np.random.randint(1,l+1)
        stop = start+PARAMS['U']-1
        #probleme si U est plus grand qu'un gene ! A corriger ! peut etre en amont lors de l'importation...
        while ((sum((l-data['TSS']['TSS_pos']+start)%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos']))) | 
            (sum((l-stop+data['TTS']['TTS_pos'])%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos']))) |
            (sum(((l-start+data['TSS']['TSS_pos'])%l>=0) & ((l-start+data['TTS']['TTS_pos'])%l<=PARAMS['U']))) |
            (sum((l-start+data['Prot']['prot_pos'])%l<PARAMS['U']))) : #Prot ok)))  
            start=np.random.randint(1,l+1)
            stop = start+PARAMS['U']-1

        #deletion 
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']-= (data['TTS']['TTS_pos']>stop)*PARAMS['U']
        updated_data['TSS']['TSS_pos']-= (data['TSS']['TSS_pos']>stop)*PARAMS['U']
        updated_data['Prot']['prot_pos'] -= (data['Prot']['prot_pos']>stop)*PARAMS['U']
        updated_data['GFF']['seq']['1'] -= (data['GFF']['seq']['1']>stop)*PARAMS['U']
        updated_data['GFF']['seq'][str(l)] -= (data['GFF']['seq']['1']>stop)*PARAMS['U'] #changer si syntaxes bizarres ?
        updated_data['GFF']['seq'].columns.values[4] = str(l-PARAMS['U']) 
        updated_data['GFF']['seq_length']=l-PARAMS['U']
        return(updated_data)

    def U_insertion(self) : 
        data = self.data
        l = data['GFF']['seq_length']
        #localisation
        start=np.random.randint(1,l+1)
        while ((sum((l-data['TSS']['TSS_pos']+start)%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos']))) 
            | (sum((l-start+data['TTS']['TTS_pos'])%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos'])))) : 
            start=np.random.randint(1,l+1)
        #insertion 
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']+= (data['TTS']['TTS_pos']>start)*PARAMS['U']
        updated_data['TSS']['TSS_pos']+= (data['TSS']['TSS_pos']>start)*PARAMS['U']
        updated_data['Prot']['prot_pos'] += (data['Prot']['prot_pos']>start)*PARAMS['U']
        updated_data['GFF']['seq']['1'] += (data['GFF']['seq']['1']>start)*PARAMS['U']
        updated_data['GFF']['seq'][str(l)] += (data['GFF']['seq']['1']>start)*PARAMS['U'] #changer si syntaxes bizarres ?
        updated_data['GFF']['seq'].columns.values[4] = str(l+PARAMS['U']) 
        updated_data['GFF']['seq_length']=l+PARAMS['U']
        return(updated_data)
