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

def start_evol_simulation(config_path) :
    
    # create plasmid with config
    config = utils.read_config(config_path)
    
    p1 = Plasmid(config)

    p1.mutate(0,0)
    p1.save()
   
    #lancer la simulation : mutation, (N=1 : création d'un fichier next_params.ini et des fichiers 
    #TTS TSS GFF Prot ... _next correspondant)
    #N!=1 : actualisation des fichiers de donné.
    # le calcul de la fitness sur params_next, sa comparaison,  au sein d'un individu. 
    #Note : faire le tri dans les outputs enregistrés de simulation.py
    #actualiser ensuite les individus (si un seul individu inclure Metropolis dans l'individu.)
    #si plusieurs individus, reprod sexuelle impossible... duplication simple ?
    #differencier les logs si N = 1 ou N > 1 
    #ptit plot de l'évolution de la fitness max qui fait plaiz'
    
    """  TODO
    
        - implement simulation recover/resume
        - implement simulation repetition
        - implement automated graph generation

    """

class Plasmid:

    def __init__(self, config, ID = 1):
        
        self.config = config
        self.CONFIG = utils.parse_config(self.config)
        
        # Some information
        
        print('Plasmid     | %s'%self['PLASMID_PATH'])
        print('Environment | %s'%self['TARGET_PATH'])
        print('Parameters  | %s'%self['CONFIG_NAME'])
        
        # creating required attributes
        
        self.ID = ID
        self.rep = -1
        self.time = -1
        self.do_my = {'Ins':self.U_insertion, 
                      'Del':self.U_deletion, 
                      'Inv':self.U_inversion } #TODO : add inversion
        
        self['PROBS'] /= np.sum(self['PROBS'])
        
        os.system('mkdir -p ' + self['WPATH']) # TODO: if doesn't exists
        
        self.target = pd.read_csv(self['TARGET_PATH'],
                                  sep = '\t',
                                  header = None, 
                                  names = ['gene_id','expression'])
            

        # setting history attributes and buildig DATA
        self.history = {'fitness':[],
                        'time':[],
                        'repetition':[],
                        'event':[],
                        'kept':[],
                        'plasmid':[],
                        'plasmid_size':[],
                        'gene_ratio':[],
                        'up_down_ratio':[],
                        'mean_space':[]
                        }
       
        if self['STEPS_DONE'] == 0 and self['REPS_DONE'] == 0 :
       
            print('Simulation  | NEW')    
       
            self.restart_plasmid()
            
        
        else :
            
            print('Simulation  | RESUME(T=%d, R=%s)'%(self['STEPS_DONE'], self['REPS_DONE']))
            
            # TODO : setup these from history
            
            self.history['fitness'] = []
            self.history['time'] = []
            self.history['event'] = []
            self.history['kept'] = []
            self.history['plasmid'] = []
            self.history['repetition'] = []

            self.history['plasmid_size'] = []
            self.history['gene_ratio'] = []
            self.history['up_down_ratio'] = []
            self.history['mean_space'] = []

            self.update_plasmid_description()
            
            self.time = self.init_time
            self.repetition = self.init_reps
        
        
        self.check_config()
        
        #sortie : time event fitness TODO : rajouter metadonnées (mean distance, etc)
        #TODO :  Charles to Baptiste : au passage, comme tu auras fais du calcul de distance entre gène, 
        # est-ce que tu
        # peux aussi calculer la taille max d'un gène et modifier le 1000 du assert suivant stp ?

    def __getitem__(self, key) :
        
        return self.CONFIG[key]
        
    def __setitem__(self, key, value) :
        
        self.CONFIG[key] = value

    def __contains__(self, key) :
        
        return key in self.CONFIG

    def mutate(self, sim = 1, rep = 1) : 

        sim = self.CONFIG['SIM_TIME'] if sim <= 0 else sim
        rep = self['N_REPS'] if rep <= 0 else rep
    
        for repetition in range(self.rep, rep) :
            
            if repetition > 0 : self.restart_plasmid()
            
            for simulation in range(self.time, sim + 1) :
                
                print('\n%sREP %d%s'%('-'*25, repetition, '-'*25))
            
                #MUTATION
                choice = np.random.choice(['Ins','Del','Inv'], p = self['PROBS']) 
                
                print('T = %d\n\tOperation:\t%s'%(simulation, choice))
            
                apply_mut = self.do_my[choice]
                updated_data = apply_mut(self.data)

                #SELECTION
                next_fitness = self.get_fitness(updated_data)
                keep_new = self.keep_mutated(next_fitness)
                
                if keep_new :
                    
                    print('\tPlasmid kept')
                    
                    self.data = updated_data
                    self.fitness = next_fitness
                    
                    utils.save_data(updated_data, self.CONFIG)
                
                # Statistics
                plasmid_size = self.data['GFF']['seq_length']
                
                gene_ratio = np.sum(np.abs(self.data['TTS']['TTS_pos'] - \
                                    self.data['TSS']['TSS_pos']))/plasmid_size
                
                UD_ratio = np.sum(self.data['TSS']['TUorient'] == '+') / \
                           np.sum(self.data['TSS']['TUorient'] == '-') if \
                           np.sum(self.data['TSS']['TUorient'] == '-') != 0 else \
                           -1
                
                starts = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1)
                stops  = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1)
                mean_space = np.mean(np.abs(stops-starts))

                # History
                self.history['time'].append(simulation)
                self.history['repetition'].append(repetition)
                self.history['event'].append(choice)
                self.history['fitness'].append(self.fitness)
                self.history['kept'].append(keep_new)

                self.history['plasmid_size'].append(plasmid_size)
                self.history['gene_ratio'].append(gene_ratio)           
                self.history['up_down_ratio'].append(UD_ratio)                
                self.history['mean_space'].append(mean_space)
                
                ### RAJOUTER UNE DESCRIPTION COMPLETE DU GENOME.
                
                self.update_plasmid_description()
                
                self.time += 1
        
    def get_fitness(self, data) :
        
        proportions = simulation.start_transcribing_2(self.config, data) 
        proportions = proportions/sum(proportions)
        
        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        
        diff = -abs(proportions-self.target['expression'].values)/len(proportions)
        
        return(round(np.sum(diff),5))

    def keep_mutated(self,next_fitness):
        print('\tFitness:\t%f -> %f'%(round(self.fitness,4), round(next_fitness,4)))
        if next_fitness>self.fitness:
            return(True)
        else : 
            alpha = math.exp(-self['ALPHA_C']*abs(self.fitness-next_fitness))
            print('\tAlpha:    \t%f'%alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))
    
    def normalize_plasmid(self):
    
        offset = self.data['TSS']['TSS_pos'][0] - 1

        self.data['GFF']['seq']['start'] -= offset
        self.data['GFF']['seq']['end'] -= offset
        self.data['Prot']['prot_pos'] = (self.data['Prot']['prot_pos']-offset)%self.data['GFF']['seq_length']

        self.data['TTS']['TTS_pos'] -= offset
        self.data['TSS']['TSS_pos'] -= offset

        return
    
    def check_config(self):

        min_gene_size = np.min(self.data['TTS']['TTS_pos']-self.data['TSS']['TSS_pos'])
        
        assert self['U'] < min_gene_size
        assert np.all(self['PROBS'] > 0)
        assert np.sum(self['PROBS']) == 1
        assert np.all(self.data['TSS']['TSS_pos'][0] > 0)
        assert np.all(self.data['TTS']['TTS_pos'][0] > 0)
        assert self.data['TSS']['TSS_pos'][0] == 1
        
        return
    
    def restart_plasmid(self) :
        
        self.rep += 1
        self.time = 0
        
        # Reading original data
        self.data = utils.read_plasmid_data(self.config)
        utils.save_data(self.data, self.CONFIG)
        self.normalize_plasmid()
        
        # Computing first line of history
        self.fitness = self.get_fitness(self.data)

        plasmid_size = self.data['GFF']['seq_length']

        gene_ratio = np.sum(np.abs(self.data['TTS']['TTS_pos'] - \
                            self.data['TSS']['TSS_pos']))/plasmid_size
            
        UD_ratio = np.sum(self.data['TSS']['TUorient'] == '+') / \
                   np.sum(self.data['TSS']['TUorient'] == '-') if \
                   np.sum(self.data['TSS']['TUorient'] == '-') != 0 else \
                   -1
        
        starts = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1)
        stops  = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1)
        mean_space = np.mean(np.abs(stops-starts))

        self.history['fitness'].append(self.fitness)
        self.history['time'].append(0)
        self.history['repetition'].append(self.rep)
        self.history['event'].append('Beg')
        self.history['kept'].append(True)
        #self.history['plasmid']
        
        self.history['plasmid_size'].append(plasmid_size)
        self.history['gene_ratio'].append(gene_ratio)
        self.history['up_down_ratio'].append(UD_ratio)
        self.history['mean_space'].append(mean_space)
        
        self.update_plasmid_description()
        
        self.time += 1
        
        return
    
    def update_plasmid_description(self):
        
        #TUindex	TUorient	TSS_pos	TSS_strength
        #TUindex	TUorient	TTS_pos	TTS_proba_off
        #prot_name	prot_pos
        #time    type    position length    strand
        
        tim = self.time
        typ = 'G'
        rep = self.rep
        
        for index in range(len(self.data['TSS'])) :
            
            pos = self.data['TSS']['TSS_pos'][index]
            ori = self.data['TSS']['TUorient'][index]
            lth = np.abs(self.data['TTS']['TTS_pos'][index] - self.data['TSS']['TSS_pos'][index])
            
            
            ori = (1 if ori == '+' else -1)
            
            self.history['plasmid'] += ['%d\t%d\t%s\t%d\t%d\t%d'%(rep, tim, typ, pos, lth, ori)]
        
        typ = 'P'
        ori = 0
        lth = 1
        
        for index in range(len(self.data['Prot'])) :
            
            pos = self.data['Prot']['prot_pos'][index]
            
            self.history['plasmid'] += ['%d\t%d\t%s\t%d\t%d\t%d'%(rep, tim, typ, pos, lth, ori)]
    
    def save_plasmid_description(self) :
        
        # Build string to save
        res = 'repetition\ttime\ttype\tlocation\tlength\tstrand\n'
        res = res + '\n'.join(self.history['plasmid'])
        
        # saving location
        location = self['PLASMID_SAVE_PATH']
        
        # Save string
        with open(location, 'w') as f :
            
            f.write(res)

        print('Plasmid configuration saved to\t%s'%location)

    def save_history(self) :
    
        # Build history data frame
        
        subdict = { key : self.history[key] for key in self.history if key != 'plasmid' }
        
        df = pd.DataFrame(data = subdict)
    
        # saving location
        location = self['HISTORY_SAVE_PATH']
        
        # Save to CSV
        df.to_csv(path_or_buf = location, index = False, sep = '\t')
    
        print('Fitness history saved to\t%s'%location)

    def U_inversion(self, data):
        
        updated_data = copy.deepcopy(data)
        
        l = data['GFF']['seq_length']     
        
        a = np.random.randint(0,l)
        b = np.random.randint(0,l)
        
        (a, b) = (b, a) if  a > b else  (a, b) # reversing if necessary

        start = updated_data['GFF']['seq'].apply(lambda x : min(x['start'],x['end']),axis=1)
        end = updated_data['GFF']['seq'].apply(lambda x : max(x['start'],x['end']),axis=1)
 
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

    def U_deletion(self, data) :
        
        updated_data = copy.deepcopy(data)
 
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
        updated_data = copy.deepcopy(data)
        updated_data['TTS']['TTS_pos']-= (data['TTS']['TTS_pos']>stop)*self['U']
        updated_data['TSS']['TSS_pos']-= (data['TSS']['TSS_pos']>stop)*self['U']
        updated_data['Prot']['prot_pos'] -= (data['Prot']['prot_pos']>stop)*self['U']
        updated_data['GFF']['seq']['start'] -= (data['GFF']['seq']['start']>stop)*self['U']
        updated_data['GFF']['seq']['end'] -= (data['GFF']['seq']['start']>stop)*self['U'] #changer si syntaxes bizarres ?
        
        updated_data['GFF']['seq_length'] = l-self['U']
        
        return(updated_data)

    def U_insertion(self,data) : 
        
        updated_data = copy.deepcopy(data)

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
        
        updated_data['GFF']['seq_length'] += self['U']
        
        return(updated_data)
    
    def save(self) :
        
        self.save_history()
        self.save_plasmid_description()
        
        return

