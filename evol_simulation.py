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
    assert PARAMS['U'] < 1000 #max gene
    PARAMS['SIM_TIME'] = 50
    PARAMS['POP_SIZE'] = 1 # on ne gere que ce cas.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  np.array([0,0,1.0, 0.0],dtype=float)
    PARAMS["alpha_c"] = 700
    PARAMS['COMPUTE_FITNESS'] = True

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
        self.do_my = {  "Ins":self.U_insertion, 
                        "Del": self.U_deletion, 
                        "Inv":self.U_inversion,
                        "None": lambda x : x} #TODO : add inversion
        self.data = utils.import_data_from_params_seq_file(PARAMS["path_params_seq"]+"params_seq.ini")
        utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_1"])        
        self.fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini", PARAMS["w_path_1"])
        
        #HISTORY
        self.hist_fitness = [self.fitness]
        self.hist_timestamp = [self.time]
        self.hist_event = ["Beg"]
        self.hist_kept = [False]
        self.genes = pd.Series(["g1","g2","g3"])
        #sortie : time event fitness TODO : rajouter metadonnées (mean distance, etc)
        #TODO :  Charles to Baptiste : au passage, comme tu auras fais du calcul de distance entre gène, 
        # est-ce que tu
        # peux aussi calculer la taille max d'un gène et modifier le 1000 du assert suivant stp ?


    def mutate(self) : 

        self.time += 1

        #MUTATION
        choice = np.random.choice(["Ins","Del","Inv","None"],p = PARAMS["probs"]) 
        
        print("/// T = %d\n\tchoice : %s"%(self.time, choice))
    
        apply_mut = self.do_my[choice]
        updated_data = apply_mut(self.data)
        utils.save_data_to_path(updated_data,PARAMS["w_path_1"])

        #SELECTION
        
        next_fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini",PARAMS["w_path_1"])
        
        keep_new = self.keep_mutated(next_fitness)
        
        if keep_new :
            print("\tmutate !")
            #UPDATE BEST
            self.data = copy.deepcopy(updated_data)
            self.fitness = next_fitness
            #UPDATE HISTORY
        
        self.hist_timestamp.append(self.time)
        self.hist_event.append(choice)
        self.hist_fitness.append(self.fitness)
        self.hist_kept.append(keep_new)
            ### RAJOUTER UNE DESCRIPTION COMPLETE DU GENOME.
        
    
    def get_fitness(self,params_file, w_path):
        proportions = simulation.start_transcribing(params_file, w_path)[5]  
        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        proportions = proportions/sum(proportions)
        diff = -abs(proportions-PARAMS['perfection']["expression"].values)/len(proportions)
        return(round(np.sum(diff),5))
        
    def keep_mutated(self,next_fitness):
        print('\t' + str(round(self.fitness,4))+" versus new : "+str(round(next_fitness,4)))
        if next_fitness>self.fitness:
            return(True)
        else : 
            alpha = math.exp(-PARAMS["alpha_c"]*abs(self.fitness-next_fitness))
            print('\t%f'%alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))
    def get_plasmid_description(self):
        #idee : 
        #self.gene = {"g1":[],"g2":[]...}
        # for g in self.genes.keys()
        #   self.gene[g].append([les infos.])
        return(0)

    def U_inversion(self, data):
        
        updated_data = copy.deepcopy(data)
    
        names=["seqid", "source", "type","start","end","score","strand","phase","attributes"]
        old_names = updated_data['GFF']['seq'].columns
        updated_data['GFF']['seq'].columns = names
        data['GFF']['seq'].columns = names
        
        l = data['GFF']['seq_length']     
        
        a = np.random.randint(0,l)
        b = np.random.randint(0,l)
        
        (a, b) = (b, a) if  a > b else  (a, b) # reversing if necessary
        
        pos_plus = data['GFF']['seq']['strand'] == '+'
        pos_moins = data['GFF']['seq']['strand'] == '-'
        
        start = copy.deepcopy(data['GFF']['seq']['start'])
        stop = copy.deepcopy(data['GFF']['seq']['end'])
        
        start.loc[pos_moins] = data['GFF']['seq']['end'].loc[pos_moins]
        stop.loc[pos_moins] = data['GFF']['seq']['start'].loc[pos_moins]
        
        while   np.any(
                        np.logical_or(
                                np.logical_and(a >= start, a <= stop),
                                np.logical_and(b >= start, b <= stop))) \
                or \
                np.any(np.logical_or(a == data['Prot']['prot_pos'], b == data['Prot']['prot_pos'])) :
        
                #print('rejected cuts %d / %d'%(a,b))
        
                a = np.random.randint(0,l)
                b = np.random.randint(0,l)
        
                (a, b) = (b, a) if  a > b else  (a, b)
                
        assert(b > a)
        
        # Genes inversions
        
        genes_in_the_middle = np.logical_and(data['GFF']['seq']['start'] > a, data['GFF']['seq']['start'] < b)
        genes_oriented_forward = np.logical_and(data['GFF']['seq']['strand'] == '+', genes_in_the_middle)
        genes_oriented_backward = np.logical_and(np.logical_not(data['GFF']['seq']['strand'] == '+'), genes_in_the_middle)
        

        new_genes_start = a + (b - data['GFF']['seq']['start'][genes_in_the_middle]) #stop becomes start
        new_genes_stop  = a + (b - data['GFF']['seq']['end'][genes_in_the_middle])
        
        # Proteins inversions
        
        proteins_in_the_middle = np.logical_and( 
                                        data['Prot']['prot_pos'] > a,
                                        data['Prot']['prot_pos'] < b )
        
        
        new_proteins_positions = a + (b - data['Prot']['prot_pos'][proteins_in_the_middle])
        
        # Update data
        
        updated_data['Prot'].loc[proteins_in_the_middle, 'prot_pos'] = new_proteins_positions
        
        updated_data['GFF']['seq'].loc[genes_in_the_middle, 'start'] = new_genes_start
        updated_data['GFF']['seq'].loc[genes_in_the_middle, 'end'] = new_genes_stop
        
        updated_data['TSS'].loc[genes_in_the_middle, 'TSS_pos'] = new_genes_start
        updated_data['TTS'].loc[genes_in_the_middle,'TTS_pos'] = new_genes_stop
        
        
        updated_data['GFF']['seq'].loc[genes_oriented_forward, 'strand'] = '-'
        updated_data['GFF']['seq'].loc[genes_oriented_backward, 'strand'] = '+'
        
        updated_data['TSS'].loc[genes_oriented_forward, 'TUorient'] = '-'
        updated_data['TSS'].loc[genes_oriented_backward, 'TUorient'] = '+'
        
        updated_data['TTS'].loc[genes_oriented_forward, 'TUorient'] = '-'
        updated_data['TTS'].loc[genes_oriented_backward, 'TUorient'] = '+'
        
        # Sorting dataframes
        
        #updated_data['GFF']['seq'].sort_values(by='start', inplace=True)
        #updated_data['TSS'].sort_values(by='TSS_pos', inplace=True)
        #updated_data['TTS'].sort_values(by='TTS_pos', inplace=True)
        #updated_data['Prot'].sort_values(by='prot_pos', inplace=True)
    
        updated_data['GFF']['seq'].columns = old_names
        data['GFF']['seq'].columns = old_names
        
        # Debug section
     
        print(data['GFF']['seq'])
        print(data['Prot'])
        #print(genes_in_the_middle)
        #print(new_genes_start)
        #print(new_genes_stop)
        
        #print('start: %d / stop: %d'%(a,b))
        #print(data['Prot'])
        #print(proteins_in_the_middle)
        #print(data['Prot']['prot_pos'][proteins_in_the_middle])
        #print(new_proteins_positions)
        
        #print('start: %d / stop: %d'%(a,b))
        #print(data['GFF']['seq'])
        
        #print(updated_data['GFF']['seq'])
        #print(updated_data['Prot'])
        
        return(updated_data)
        
    def U_deletion(self,data) :
        
        updated_data = copy.deepcopy(data)
        
        names = ["seqid", "source", "type","start","end","score","strand","phase","attributes"]
        old_names = data['GFF']['seq'].columns
        
        updated_data['GFF']['seq'].columns = names
        data['GFF']['seq'].columns = names
        
        l = data['GFF']['seq_length']     
        
        #localisation
        start=np.random.randint(0,l-PARAMS['U'])
        stop = start+PARAMS['U']-1 #OKKKK : hyp que premiere base est premier gene.
        #sinon relancer
        #probleme si U est plus grand qu'un gene ! A corriger ! peut etre en amont lors de l'importation...
        
        while (sum( ((data['TSS']['TSS_pos'] - start) >= PARAMS['U']))): #Prot ok)))  
            start=np.random.randint(0,l-PARAMS['U'])
            stop = start+PARAMS['U']-1 

        #deletion 

        #deletion 
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']-= (data['TTS']['TTS_pos']>stop)*PARAMS['U']
        updated_data['TSS']['TSS_pos']-= (data['TSS']['TSS_pos']>stop)*PARAMS['U']
        updated_data['Prot']['prot_pos'] -= (data['Prot']['prot_pos']>stop)*PARAMS['U']
        updated_data['GFF']['seq']['start'] -= (data['GFF']['seq']['start']>stop)*PARAMS['U']
        updated_data['GFF']['seq']['end'] -= (data['GFF']['seq']['start']>stop)*PARAMS['U'] #changer si syntaxes bizarres ?
        
        new_header = old_names.values
        new_header[4] = str(l-PARAMS['U'])
        
        updated_data['GFF']['seq_length'] = l-PARAMS['U']
        updated_data['GFF']['seq'].columns = new_header
        
        data['GFF']['seq'].columns = new_header
        
        return(updated_data)

    def U_insertion(self,data) : 
        
        updated_data = copy.deepcopy(data)
        
        names=["seqid", "source", "type","start","end","score","strand","phase","attributes"]
        old_names = data['GFF']['seq'].columns
        
        updated_data['GFF']['seq'].columns = names
        data['GFF']['seq'].columns = names
        
        data = self.data
        l = data['GFF']['seq_length']
        #localisation
        start=np.random.randint(1,l+1)
        while ((sum((l-data['TSS']['TSS_pos']+start)%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos']))) 
            | (sum((l-start+data['TTS']['TTS_pos'])%l<abs(data['TTS']['TTS_pos']-data['TSS']['TSS_pos'])))) : 
            start=np.random.randint(1,l+1)
        #insertion 

        updated_data['TTS']['TTS_pos']+= (data['TTS']['TTS_pos']>start)*PARAMS['U']
        updated_data['TSS']['TSS_pos']+= (data['TSS']['TSS_pos']>start)*PARAMS['U']
        updated_data['Prot']['prot_pos'] += (data['Prot']['prot_pos']>start)*PARAMS['U']
 
        updated_data['GFF']['seq']['start'] += (data['GFF']['seq']['start']>start)*PARAMS['U']
        updated_data['GFF']['seq']['end'] += (data['GFF']['seq']['start']>start)*PARAMS['U'] 

        new_header = old_names.values
        new_header[4] = str(l+PARAMS['U'])
        
        updated_data['GFF']['seq_length']+=PARAMS['U']
        
        updated_data['GFF']['seq'].columns = new_header
        data['GFF']['seq'].columns = new_header
        
        return(updated_data)
    
    def save(self) :
        
        # Save history
        
        self.history = pd.DataFrame(data={  'time' : self.hist_timestamp,
                                            'event' : self.hist_event, 
                                            'fitness' : self.hist_fitness,
                                            'kept' : self.hist_kept})
        
        path = PARAMS["w_path_1"]+"history.csv"
        print("Saving evolution history to : "+path)
        
        self.history.to_csv(path_or_buf = path, index=False, sep='\t',
                            columns=['time', 'event', 'kept', 'fitness'])
    
        # Save
    
        utils.save_data_to_path(self.data, PARAMS["w_path_1"] + '_final')

