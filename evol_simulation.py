# coding: utf-8
import pandas as pd
import numpy as np
import copy
import os
import simulation
import utils
PARAMS = {}

"""
TODO :
- parser au debut de start_evol_simulation (lister au passage paramères dans readme.md)
- Problème avec le parser dans start_transcribing lorsqu'on est dans le working dir
-> fonction de utils ne trouve pas de section ; probablement un problème de path qui n'est pas
absolu. Essayer de résoudre !
- Qd n!=1 :
#NONcréer une population et la stocker
        #chaque individu est crée avec le fichier de params.ini comme argument, #PAS ENCORE car N=1

    #boucler SIM_TIME fois (en nombre de mut ?)
        #tirer les mutations possibles et les exécuter ; 
"""

def start_evol_simulation(INI_file) :

    #import model parameters from INI_fileparams_evol.init (dont le nom du fichier params.init)
    #note : utiliser syntaxe de simulation.read_config_file(path) ?
    PARAMS['U'] = 150
    PARAMS['POP_SIZE'] = 1 # on ne gere que ce cas.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  [1/3.0,1/3.0,1/3.0]
    PARAMS['perfection'] = [1/10.0]*10
    PARAMS['w_path_1'] = "tousgenesidentiques/N_1/pIns_0.33/pDel_0.33/pInv_0.33/1/" #TODO : automatiser
    PARAMS['w_path_2'] = "tousgenesidentiques/N_1/pIns_0.33/pDel_0.33/pInv_0.33/2/"
    os.system("mkdir -p " + PARAMS['w_path_1']) #à adapter et peut etre deplacer dans plasmid.
    os.system("mkdir -p " + PARAMS['w_path_2'])

    Plasmid = plasmid()
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
        self.time = 0
        self.do_my = {"Ins":self.U_insertion, "Del": self.U_deletion, "Inv":self.U_inversion} #TODO : add inversion
        self.data = utils.import_data_from_params_seq_file(PARAMS["path_params_seq"]+"params_seq.ini")
        
        # copy of files in the working paths
        utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_1"])
        utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_2"])
        
        #self.fitness = self.get_fitness(PARAMS["w_path_1"]+"params.seq", PARAMS["w_path_1"])
        self.fitness = 500 #TODO : change
        self.hist_fitness = [self.fitness]
        self.hist_event = []
        #TODO : réfléchir à l'historique dans le cas ou on n'a pas forcément un evenement qui arrive !
        # possibilité pour event : créer une ligne text de la forme
        # Inv 1546
        # (seule)    

    def mutate(self) : 

        choice = np.random.choice(["Ins","Del","Inv"],p = PARAMS["probs"]) #pour linstant systématique. Changer ?

        #apply_mut = self.do_my[choice]
        apply_mut = self.U_insertion #for testing purpose


        print("The event is : " +choice)
        #print(self.data)
        updated_data = apply_mut(self.data)
        #print(updated_data)
        utils.save_data_to_path(updated_data,PARAMS["w_path_2"])

        #next_fitness = self.get_fitness(PARAMS["w_path_2"]+"params.seq",PARAMS["w_path_2"])
        next_fitness = 1000 #TODO : change
        if self.keep_mutated(next_fitness) :

            self.fitness = next_fitness
            self.hist_fitness.append(self.fitness)
            self.hist_event.append(choice)
            self.data = copy.deepcopy(updated_data)  #Deep-copy obligatoire ?
            
            #chosir enfin d'alterner les working space ou non...
    
    def get_fitness(self,params_file, w_path):
        output = simulation.start_transcribing(params_file, w_path)
        nb_transcribed = np.array([10.0,20.0],dtype=float) 
        #TODO : trouver lequel element c'est dans output et si les identifiants correspondent.
        proportions = nb_transcribed/sum(nb_transcribed)
        return(-abs(proportions-PARAMS["perfection"])/len(nb_transcribed))
        
    def keep_mutated(self,next_fitness):
        if next_fitness>self.fitness:
            return(True)
        else : 
            e = self.fitness-next_fitness
            return(np.random.choice([True,False],p = [0.1,0.9])) ## TODO : write formula for p !!

    #TODO
    def U_inversion(self,data):
        
        return(copy.deepcopy(data))
        
    def U_deletion(self,data) :
        
        l = data['GFF']['seq_length']
        
        print(PARAMS['U'])
        
        #localisation
        start=np.random.randint(1,l+1)
        stop = start+PARAMS['U']-1
        #probleme si U est plus grand qu'un gene ! A corriger ! peut etre en amont lors de l'importation...
        
    
        # Why (x + l) %l  ??
        #while ( (sum( (l - data['TSS']['TSS_pos'] + start) % l < abs( data['TTS']['TTS_pos'] - data['TSS']['TSS_pos'] ) ) )     | 
                #(sum( (l + data['TTS']['TTS_pos'] - stop)  % l < abs( data['TTS']['TTS_pos'] - data['TSS']['TSS_pos'] ) ) )      |
                #(sum( ((l - start + data['TSS']['TSS_pos']) % l >= 0) & ((l - start+data['TTS']['TTS_pos']) % l <= PARAMS['U'] ) ) ) |
                #(sum( (l - start + data['Prot']['prot_pos']) % l < PARAMS['U'] ) ) ) : #Prot ok)))  
            
            #start=np.random.randint(1,l+1)
            
            #stop = start+PARAMS['U']-1

        while sum( ( ( start > data['TSS']['TSS_pos']  ) & ( start  < data['TSS']['TSS_pos'] ) ) |
                ( ( stop  > data['TSS']['TSS_pos']  ) & ( stop   < data['TSS']['TSS_pos'] ) ) |
                ( ( data['Prot']['prot_pos'] > data['TSS']['TSS_pos']  ) & ( data['Prot']['prot_pos'] < data['TSS']['TSS_pos'] ) ) ) :
                
            
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

    def U_insertion(self,data) : 
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
    
    def save(self) :
        #to save hist_fitness, and maybe other stuff, later
        pass
    
