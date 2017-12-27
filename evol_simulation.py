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



"""
TODO :
- parser au debut de start_evol_simulation (lister au passage paramères dans readme.md)
"""

def start_evol_simulation(config_path) :
    
    # create plasmid with config with config
    CONFIG = utils.read_config(config_path)
    Plasmid = plasmid(CONFIG)

    Plasmid.mutate()
    Plasmid.save()

    # if self['path_params_seq'] == "tousgenesidentiques/":
    #     fit_std = 0.000758179842093
    # else : 
    #     fitness_base = []
    #     for i in range(10):
    #         fitness_base.append(Plasmid.get_fitness(self["w_path_1"]+"params_seq.ini",self["w_path_1"]))
    #     print(fitness_base)
    #     fitness_base = np.array(fitness_base)
    #     fit_std = np.std(fitness_base)
    #     fit_mean = np.mean(fitness_base)
    #     n_rep_1_perc = math.ceil((fit_std/(0.01*fit_mean))**2)
    #     print(fit_std, fit_mean,n_rep_1_perc)

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

    def __init__(self, config, ID = 1, history_file = 'history', plasmid_file = 'plasmid'):
        
        self.config = config
        
        # cheking config
        assert self['U'] < 1000          #max gene TODO: automatiser
        assert np.all(self['PROBS'] > 0)

        # creating required attributes
        self['PROBS'] /= np.sum(self['PROBS'])
        os.system("mkdir -p " + self['WPATH'])
        self.target = pd.read_csv(#'tousgenesidentiques/environment.dat',
                                    self['TARGET_PATH'],
                                  sep = '\t',
                                  header = None, 
                                  names = ["gene_id","expression"])
        self.data = utils.import_data_from_params_seq_file(self['PARAMS_PATH'])
        
        self.ID = ID
        self.time = 0
        self['HISTORY_PATH'] = self['WPATH'] + history_file + '.csv'
        self['PLASMID_PATH'] = self['WPATH'] + plasmid_file + '.csv'
        self.do_my = {"Ins":self.U_insertion, 
                      "Del":self.U_deletion, 
                      "Inv":self.U_inversion } #TODO : add inversion

        # normalizing plamid TODO : AUTOMATISER
        self.data['GFF']['seq']["1"]-=1000
        self.data['GFF']['seq']["30000"]-=1000
        self.data['Prot']['prot_pos']= (self.data['Prot']['prot_pos']-1000)%self.data["GFF"]["seq_length"]

        self.data['TTS']["TTS_pos"]-=1000
        self.data['TSS']["TSS_pos"]-=1000


        self.fitness = self.get_fitness(self['PARAMS_PATH'],
                                        self['WPATH'])  
        
        # setting history attributes
        self.hist_fitness = [self.fitness]
        self.hist_timestamp = [self.time]
        self.hist_event = ["Beg"]
        self.hist_kept = [False]
        self.hist_plasmid = []
        
        self.update_plasmid_description()  
        #sortie : time event fitness TODO : rajouter metadonnées (mean distance, etc)
        #TODO :  Charles to Baptiste : au passage, comme tu auras fais du calcul de distance entre gène, 
        # est-ce que tu
        # peux aussi calculer la taille max d'un gène et modifier le 1000 du assert suivant stp ?

    def __getitem__(self, key) :
        
        return self.config[key]
        
    def __setitem__(self, key, value) :
        
        self.config[key] = value

    def __contains__(self, key) :
        
        return key in self.config

    def mutate(self) : 

        self.time += 1

        #MUTATION
        choice = np.random.choice(["Ins","Del","Inv"],p = self.config['PROBS']) 
        
        print("/// T = %d\n\tchoice : %s"%(self.time, choice))
    
        apply_mut = self.do_my[choice]
        updated_data = apply_mut(self.data)

        #SELECTION
        
        next_fitness = self.get_fitness(self['PARAMS_PATH'], 
                                        self['WPATH'])
        
        keep_new = self.keep_mutated(next_fitness)
        
        if keep_new :
            print("\tmutate !")
            #UPDATE BEST
            utils.save_data_to_path(updated_data, self['WPATH'])
            self.data = copy.deepcopy(updated_data)
            self.fitness = next_fitness
            #UPDATE HISTORY
        
        self.hist_timestamp.append(self.time)
        self.hist_event.append(choice)
        self.hist_fitness.append(self.fitness)
        self.hist_kept.append(keep_new)
            ### RAJOUTER UNE DESCRIPTION COMPLETE DU GENOME.
        
        self.update_plasmid_description()
    
    def get_fitness(self, params_file, w_path) :
        
        proportions = simulation.start_transcribing_2(params_file,self.data, w_path) 
        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        proportions = proportions/sum(proportions)
        diff = -abs(proportions-self.target["expression"].values)/len(proportions)
        return(round(np.sum(diff),5))
        
    def keep_mutated(self,next_fitness):
        print('\t' + str(round(self.fitness,4))+" versus new : "+str(round(next_fitness,4)))
        if next_fitness>self.fitness:
            return(True)
        else : 
            alpha = math.exp(-self['ALPHA_C']*abs(self.fitness-next_fitness))
            print('\t%f'%alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))
    
    '''
        Warning : do not use twice in a simulation
    '''
    def update_plasmid_description(self):
        
        #TUindex	TUorient	TSS_pos	TSS_strength
        #TUindex	TUorient	TTS_pos	TTS_proba_off
        #prot_name	prot_pos
        #time    type    position length    strand
        
        tim = self.time
        typ = 'G'
        
        for index in range(len(self.data['TSS'])) :
            
            pos = self.data['TSS']['TSS_pos'][index]
            ori = self.data['TSS']['TUorient'][index]
            lth = np.abs(self.data['TTS']['TTS_pos'][index] - self.data['TSS']['TSS_pos'][index])
            
            ori = (1 if ori == '+' else -1)
            
            self.hist_plasmid += ['%d\t%s\t%d\t%d\t%d'%(tim, typ, pos, lth, ori)]
        
        typ = 'P'
        ori = 0
        lth = 1
        
        for index in range(len(self.data['Prot'])) :
            
            pos = self.data['Prot']['prot_pos'][index]
            
            self.hist_plasmid += ['%d\t%s\t%d\t%d\t%d'%(tim, typ, pos, lth, ori)]
    
    
    def save_plasmid_description(self, location) :
        
        # Build string to save
        res = 'time\ttype\tlocation\tlength\tstrand\n'
        
        res = res + '\n'.join(self.hist_plasmid)
        
        # Save string
        with open(location, 'w') as f :
            
            f.write(res)
            
        print('Plasmid configuration saved to \'%s\''%location)
        
        
    def save_history(self, location) :
    
        # Build history data frame
        self.history = pd.DataFrame(data={  'time' : self.hist_timestamp,
                                            'event' : self.hist_event, 
                                            'fitness' : self.hist_fitness,
                                            'kept' : self.hist_kept})
        
        # Save to CSV
        self.history.to_csv(path_or_buf = location, index = False, sep = '\t',
                            columns = ['time', 'event', 'kept', 'fitness'])
    
        print('Fitness history saved to %s'%location)


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

        start = updated_data["GFF"]["seq"].apply(lambda x : min(x["start"],x["end"]),axis=1)
        end = updated_data["GFF"]["seq"].apply(lambda x : max(x["start"],x["end"]),axis=1)
        print(start,end)
        while   np.any(
                        np.logical_or(
                                np.logical_and(a >= start, a<=end),
                                np.logical_and(b >= start, b<=end))) \
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
     
        #print(data['GFF']['seq'])
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
        start=np.random.randint(0,l-self['U'])
        stop = start+self['U']-1 #OKKKK : hyp que premiere base est premier gene.
        #sinon relancer
        #probleme si U est plus grand qu'un gene ! A corriger ! peut etre en amont lors de l'importation...
        
        while (sum( ((data['TSS']['TSS_pos'] - start) >= self['U']))): #Prot ok)))  
            start=np.random.randint(0,l-self['U'])
            stop = start+self['U']-1 

        #deletion 

        #deletion 
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']-= (data['TTS']['TTS_pos']>stop)*self['U']
        updated_data['TSS']['TSS_pos']-= (data['TSS']['TSS_pos']>stop)*self['U']
        updated_data['Prot']['prot_pos'] -= (data['Prot']['prot_pos']>stop)*self['U']
        updated_data['GFF']['seq']['start'] -= (data['GFF']['seq']['start']>stop)*self['U']
        updated_data['GFF']['seq']['end'] -= (data['GFF']['seq']['start']>stop)*self['U'] #changer si syntaxes bizarres ?
        
        new_header = old_names.values
        new_header[4] = str(l-self['U'])
        
        updated_data['GFF']['seq_length'] = l-self['U']
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

        updated_data['TTS']['TTS_pos']+= (data['TTS']['TTS_pos']>start)*self['U']
        updated_data['TSS']['TSS_pos']+= (data['TSS']['TSS_pos']>start)*self['U']
        updated_data['Prot']['prot_pos'] += (data['Prot']['prot_pos']>start)*self['U']
 
        updated_data['GFF']['seq']['start'] += (data['GFF']['seq']['start']>start)*self['U']
        updated_data['GFF']['seq']['end'] += (data['GFF']['seq']['start']>start)*self['U'] 

        new_header = old_names.values
        new_header[4] = str(l+self['U'])
        
        updated_data['GFF']['seq_length']+=self['U']
        
        updated_data['GFF']['seq'].columns = new_header
        data['GFF']['seq'].columns = new_header
        
        return(updated_data)
    
    def save(self) :
        
        self.save_history(self['HISTORY_PATH'])
        self.save_plasmid_description(self['PLASMID_PATH'])



