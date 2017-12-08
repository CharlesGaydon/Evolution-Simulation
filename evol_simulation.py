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
"""

def start_evol_simulation(INI_file) :

    #import model parameters from INI_fileparams_evol.init (dont le nom du fichier params.init)
    #note : utiliser syntaxe de simulation.read_config_file(path) ?
    PARAMS['U'] = 150
    assert PARAMS['U']<1000 #max gene
    PARAMS['SIM_TIME'] = 50
    PARAMS['POP_SIZE'] = 1 # on ne gere que ce cas.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  [0.0, 0.0, 1.0]

    assert(sum(PARAMS['probs'])<=1) 

    PARAMS['w_path_1'] = "tousgenesidentiques/pIns_0.33/pDel_0.33/pInv_0.33/" #TODO : automatiser

    PARAMS['perfection'] = pd.read_csv(filepath_or_buffer = PARAMS['path_params_seq']+"environment.dat",
                        sep = '\t',
                        header=None, 
                        names = ["gene_id","expression"])

    os.system("mkdir -p " + PARAMS['w_path_1']) #à adapter et peut etre deplacer dans plasmid.


    Plasmid = plasmid()
    while Plasmid.time<PARAMS['SIM_TIME']:
        Plasmid.mutate()
    Plasmid.save()

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
        #TODO 
        # copy of files in the working paths
        utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_1"])
        # utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_2"])
        
        self.fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini", PARAMS["w_path_1"])
        
        #HISTORY
        self.hist_fitness = [self.fitness]
        self.hist_timestamp = [self.time]
        self.hist_event = ["Beg"]
                        #    sortie : time event fitness


    def mutate(self) : 
        #MUTATION
        while self.time<=PARAMS["SIM_TIME"] :
            choice = np.random.choice(["Ins","Del","Inv","No_event"],p = PARAMS["probs"]) 
            self.time+=1
            if choice!="No_event":
                break
        else:
            return("Simulation time is exceeded so we get out of mutate!")
        #apply_mut = self.do_my[choice]
        choice = "Ins"
        print("t = : "+str(self.time)+", Choice = : " +choice)

        apply_mut = self.U_insertion #for testing purpose
        updated_data = apply_mut(self.data)
        utils.save_data_to_path(updated_data,PARAMS["w_path_1"])
        next_fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini",PARAMS["w_path_1"])

         #SELECTION
        if self.keep_mutated(next_fitness) :
            print("mutate !")
            #UPDATE BEST
            self.data = copy.deepcopy(updated_data)
            self.fitness = next_fitness
            #UPDATE HISTORY
            self.hist_timestamp.append(self.time)
            self.hist_event.append(choice)
            self.hist_fitness.append(self.fitness)
            ### RAJOUTER UNE DESCRIPTION COMPLETE DU GENOME.
    
    def get_fitness(self,params_file, w_path):
        proportions = simulation.start_transcribing(params_file, w_path)[5] 
        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        proportions = proportions/sum(proportions)
        diff = -abs(proportions-PARAMS['perfection']["expression"].values)/len(proportions)
        return(np.sum(diff))
        
    def keep_mutated(self,next_fitness):
        print(str(self.fitness)+" versus new : "+str(next_fitness))
        if next_fitness>self.fitness:
            return(True)
        else : 
            e = self.fitness-next_fitness
            return(np.random.choice([True,False],p = [0.2,0.8])) ## TODO : write formula for p !!

    #TODO
    def U_inversion(self,data):
        
        l = data['GFF']['seq_length'] 
        
        print
        
        b1 = np.random.randint(1,l+1)
        b2 = np.random.randint(1,l+1)
        
        while condition :
            
            b1 = np.random.randint(1,l+1)
            b2 = np.random.randint(1,l+1)    
        
        # inverser perfection au passage !
        return(copy.deepcopy(data))
        
    def U_deletion(self,data) :
        
        l = data['GFF']['seq_length']     
        #localisation
        start=np.random.randint(1,l+1-PARAMS['U'])
        stop = start+PARAMS['U']-1 #OKKKK : hyp que premiere base est premier gene.
        #sinon relancer
        #probleme si U est plus grand qu'un gene ! A corriger ! peut etre en amont lors de l'importation...
        
    
        # while ( (sum( (l - data['TSS']['TSS_pos'] + start) % l < abs( data['TTS']['TTS_pos'] - data['TSS']['TSS_pos'] ) ) )     | 
                # (sum( (l + data['TTS']['TTS_pos'] - stop)  % l < abs( data['TTS']['TTS_pos'] - data['TSS']['TSS_pos'] ) ) )      |
        while (sum( ((data['TSS']['TSS_pos'] - start) >= PARAMS['U']))): #Prot ok)))  
            start=np.random.randint(1,l+1-PARAMS['U'])
            stop = start+PARAMS['U']-1 

        # while sum( ( ( start > data['TSS']['TSS_pos']  ) & ( start  < data['TSS']['TSS_pos'] ) ) |
        #         ( ( stop  > data['TSS']['TSS_pos']  ) & ( stop   < data['TSS']['TSS_pos'] ) ) |
        #         ( ( data['Prot']['prot_pos'] > data['TSS']['TSS_pos']  ) & ( data['Prot']['prot_pos'] < data['TSS']['TSS_pos'] ) ) ) :
                
            
            start=np.random.randint(1,l+1)
            
            stop = start+PARAMS['U']-1

        #deletion 
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']-= (data['TTS']['TTS_pos']>stop)*PARAMS['U']
        updated_data['TSS']['TSS_pos']-= (data['TSS']['TSS_pos']>stop)*PARAMS['U']
        updated_data['Prot']['prot_pos'] -= (data['Prot']['prot_pos']>stop)*PARAMS['U']
        updated_data['GFF']['seq']['1'] -= (data['GFF']['seq']['1']>stop)*PARAMS['U']
        updated_data['GFF']['seq'][str(l)] -= (data['GFF']['seq']['1']>stop)*PARAMS['U'] #changer si syntaxes bizarres ?
        
        new_header = copy.deepcopy(data["GFF"]['seq'].columns.values)
        new_header[4] = str(l-PARAMS['U'])
        updated_data['GFF']['seq'].columns = new_header 
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
        #print(updated_data)
        updated_data['GFF']['seq']['1'] += (data['GFF']['seq']['1']>start)*PARAMS['U']
        updated_data['GFF']['seq'][str(l)] += (data['GFF']['seq']['1']>start)*PARAMS['U'] 

        new_header = copy.deepcopy(data["GFF"]['seq'].columns.values)
        new_header[4] = str(l+PARAMS['U'])
        updated_data['GFF']['seq'].columns = new_header #sinon même header car numpy object pas deep copié !  
        updated_data['GFF']['seq_length']+=PARAMS['U']
        return(updated_data)
    
    def save(self) :
        #to save hist_fitness, and maybe other stuff, later
        df = pd.DataFrame(np.transpose(np.vstack((self.hist_timestamp,self.hist_event,self.hist_fitness))),
            columns = ["Timestamp", "Event", "Fitness"])
        path = PARAMS["w_path_1"]+"history.csv"
        print("Saving evolution history to : "+path)
        df.to_csv(path_or_buf = path, index=False)
    
