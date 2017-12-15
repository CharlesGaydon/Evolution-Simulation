# coding: utf-8
import pandas as pd
import numpy as np
import copy
import os
import simulation
import utils
import math
'''
Authors : Charles GAYDON, Baptiste LAC
Oct 2017- February 2018
'''



PARAMS = {}

"""
TODO :
- parser au debut de start_evol_simulation (lister au passage paramères dans readme.md)
"""

def start_evol_simulation(INI_file) :

    #import model parameters from INI_fileparams_evol.init (dont le nom du fichier params.init)
    #note : utiliser syntaxe de simulation.read_config_file(path) ?
    
    ## IMPORT PARAMS
    PARAMS['U'] = 150
    PARAMS['SIM_TIME'] = 50
    PARAMS['POP_SIZE'] = 1 # on ne gere que ce cas.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  np.array([1/3.0,1/3.0,1/3.0],dtype=float)
    PARAMS["alpha_c"] = 700

    ## FURTHERMORE
    PARAMS['probs']/=sum(PARAMS['probs'])
    PARAMS['w_path_1'] = PARAMS['path_params_seq']+ utils.make_w_path(PARAMS['probs'])  
    os.system("mkdir -p " + PARAMS['w_path_1']) #à adapter et peut etre deplacer dans plasmid.
    PARAMS['perfection'] = pd.read_csv(filepath_or_buffer = PARAMS['path_params_seq']+"environment.dat",
                        sep = '\t',
                        header=None, 
                        names = ["gene_id","expression"])

    ## GENERATE PLASMID
    Plasmid = plasmid()

    # if PARAMS['path_params_seq'] == "tousgenesidentiques/":
    #     fit_std = 0.000758179842093
    # else : 
    #     fitness_base = []
    #     for i in range(10):
    #         fitness_base.append(Plasmid.get_fitness(PARAMS["w_path_1"]+"params_seq.ini",PARAMS["w_path_1"]))
    #     print(fitness_base)
    #     fitness_base = np.array(fitness_base)
    #     fit_std = np.std(fitness_base)
    #     fit_mean = np.mean(fitness_base)
    #     n_rep_1_perc = math.ceil((fit_std/(0.01*fit_mean))**2)
    #     print(fit_std, fit_mean,n_rep_1_perc)


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
        self.do_my = {"Ins":self.U_insertion, "Del": self.U_deletion, "Inv":self.U_inversion}
        self.data = utils.import_data_from_params_seq_file(PARAMS["path_params_seq"]+"params_seq.ini")
        utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_1"])        
        self.fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini", PARAMS["w_path_1"])
        
        #HISTORY
        self.hist_fitness = [self.fitness]
        self.hist_timestamp = [self.time]
        self.hist_event = ["Beg"]
        self.genes = pd.Series(["g1","g2","g3"])
        #sortie : time event fitness TODO : rajouter metadonnées (mean distance, etc)
        #TODO :  Charles to Baptiste : au passage, comme tu auras fais du calcul de distance entre gène, 
        # est-ce que tu
        # peux aussi calculer la taille max d'un gène et modifier le 1000 du assert suivant stp ?

        assert PARAMS["U"] < 1000, "Sequence unit 'U' for InDel should be smaller than the smallest gene. "

    def mutate(self) : 
        #MUTATION
        choice = np.random.choice(["Ins","Del","Inv"],p = PARAMS["probs"]) 
        self.time+=1
        apply_mut = self.do_my[choice] 
        print("t = : "+str(self.time)+", Choice = : " +choice)
        apply_mut = self.U_deletion #for testing purpose
        
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
        return(round(np.sum(diff),5))
        
    def keep_mutated(self,next_fitness):
        print(str(round(self.fitness,4))+" versus new : "+str(round(next_fitness,4)))
        if next_fitness>self.fitness:
            return(True)
        else : 
            alpha = math.exp(-PARAMS["alpha_c"]*abs(self.fitness-next_fitness))
            print(alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))
    def get_plasmid_description(self):
        #idee : 
        #self.gene = {"g1":[],"g2":[]...}
        # for g in self.genes.keys()
        #   self.gene[g].append([les infos.])
        return(0)

    #TODO
    def U_inversion(self,data):
        
        l = data['GFF']['seq_length'] 
        local_max_iter = 0        
        b1 = np.random.randint(1,l+1)
        b2 = np.random.randint(1,l+1)
        
        while condition :
            b1 = np.random.randint(1,l+1)
            b2 = np.random.randint(1,l+1)    
            local_max_iter+=1
            if local_max_iter>50: 
                sys.exit("A new inversion cannot be done on the genome.")
        
        # inverser perfection au passage !
        return(copy.deepcopy(data))
        
    def U_deletion(self,data) :
        
        l = data['GFF']['seq_length']     
        #localisation
        start=np.random.randint(0,l-PARAMS['U'])
        stop = start+PARAMS['U']
        local_max_iter = 0
        while sum(abs(data['TSS']['TSS_pos'].values - start) <= PARAMS['U']\
                | (abs(data['TTS']['TTS_pos'].values - stop) <= PARAMS['U'])
                | (abs(data['TSS']['TSS_pos'].values - start) + abs(data['TTS']['TTS_pos'].values - stop)\
                    == data['TTS']['TTS_pos'].values-data['TSS']['TSS_pos'].values)\
                |  (abs(data['Prot']['prot_pos'].values - start) + abs(data['Prot']['prot_pos'].values - stop)\
                    == PARAMS['U'] )): 
            start=np.random.randint(1,l+1-PARAMS['U'])
            stop = start+PARAMS['U']
            local_max_iter+=1
            if local_max_iter>50:
                sys.exit("A new deletion cannot be done on the genome.")

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
        local_max_iter = 0
        #localisation
        start=np.random.randint(1,l)
        # while sum(((data['TSS']['TSS_pos']<start)*1) *  ((data['TTS']['TTS_pos']>start)*1)):
        while sum((data['TSS']['TSS_pos'].values<start) &  (data['TTS']['TTS_pos'].values>start)):
            start=np.random.randint(1,l)
            local_max_iter+=1
            if local_max_iter>50:
                sys.exit("A new insertion cannot be done on the genome.")
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
    
