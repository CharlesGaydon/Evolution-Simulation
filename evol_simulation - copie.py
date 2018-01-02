# coding: utf-8
import pandas as pd
import numpy as np
import copy
import os
import simulation
import utils
import math
<<<<<<< HEAD
import time
=======

LOG = False

>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
'''
Authors : Charles GAYDON, Baptiste LAC
Oct 2017- February 2018
'''

def start_evol_simulation(config_path) :
    
    # create plasmid with config
    config = utils.read_config(config_path)
    
    p1 = Plasmid(config)

    p1.mutate(0,0)
    
    #lancer la simulation : mutation, (N=1 : création d'un fichier next_params.ini et des fichiers 
    #TTS TSS GFF Prot ... _next correspondant)
    #N!=1 : actualisation des fichiers de donné.
    # le calcul de la fitness sur params_next, sa comparaison,  au sein d'un individu. 
    #Note : faire le tri dans les outputs enregistrés de simulation.py
    #actualiser ensuite les individus (si un seul individu inclure Metropolis dans l'individu.)
    #si plusieurs individus, reprod sexuelle impossible... duplication simple ?
    #differencier les logs si N = 1 ou N > 1 
    
class Plasmid:

    def __init__(self, config):
        
        # creating required attributes
        self.rep = -1
        self.time = -1
        self.do_my = {'Ins':self.U_insertion, 
                      'Del':self.U_deletion, 
                      'Inv':self.U_inversion } #TODO : add inversion
        
        self.config = config
        self.CONFIG = utils.parse_config(self.config)
    
<<<<<<< HEAD
    ## IMPORT PARAMS
    PARAMS['U'] = 150
    assert PARAMS['U']<1000 #max gene
    PARAMS['SIM_TIME'] = 200
    PARAMS['POP_SIZE'] = 1 # on ne gere que ce cas.
    PARAMS['path_params_seq'] = "tousgenesidentiques/"
    PARAMS['probs'] =  np.array([1/2,1/2,0.0],dtype=float)
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
        print(Plasmid.data)
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
        self.data = utils.import_data_from_params_seq_file(PARAMS["path_params_seq"]+"params_seq.ini") #MODIFIED ! 
        #utils.copy_to_working_path(PARAMS['path_params_seq'], PARAMS["w_path_1"])
        # print(self.data)

        # AUTOMATISER
        self.data['GFF']['seq']["1"]-=1000
        self.data['GFF']['seq']["30000"]-=1000
        self.data['Prot']['prot_pos']= (self.data['Prot']['prot_pos']-1000)%self.data["GFF"]["seq_length"]

        self.data['TTS']["TTS_pos"]-=1000
        self.data['TSS']["TSS_pos"]-=1000
        # print("Params file : " + PARAMS["path_params_seq"]+"params_seq.ini")
        self.fitness = self.get_fitness(PARAMS["path_params_seq"]+"params_seq.ini")
        
        #HISTORY
        self.hist_fitness = [self.fitness]
        self.hist_timestamp = [self.time]
        self.hist_event = ["Beg"]
        self.hist_kept = [False]
        self.genes = self.initialize_genes_description()
        self.genes_after_inversion = None

        print(self.genes)
=======
        self['PROBS'] /= np.sum(self['PROBS'])
        
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
       
        # new simulation starts
        if self['SIMS_DONE'] == 0 and self['REPS_DONE'] == 0 :
       
            print('Simulation  | NEW (S: %d/%d R: %d/%d)'%(self['SIMS_DONE'],
                                                             self['N_SIM'],
                                                             self['REPS_DONE'],
                                                             self['N_REP']))
       
            os.system('mkdir -p ' + self['WPATH']) # TODO: if doesn't exists       
       
            self.restart_plasmid()
            
            self.save()
        
        # simulation is already finished
        elif self['SIMS_DONE'] >= self['N_SIM'] and self['REPS_DONE'] >= self['N_REP'] :
            
            print('Simulation  | ENDED')
        
            self.resume_plasmid()
        
        # resuming simulation
        else :
            
            print('Simulation  | RESUME (S: %d/%d R: %d/%d)'%(self['SIMS_DONE'],
                                                             self['N_SIM'],
                                                             self['REPS_DONE'],
                                                             self['N_REP']))
            
            # Update local path
            
            self.resume_plasmid()
        
        self.check_config()
        
        print('Plasmid     | %s'%self['PLASMID_PATH'])
        print('Environment | %s'%self['TARGET_PATH'])
        print('Parameters  | %s'%self['CONFIG_NAME'])
        
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
        #sortie : time event fitness TODO : rajouter metadonnées (mean distance, etc)
        #TODO :  Charles to Baptiste : au passage, comme tu auras fais du calcul de distance entre gène, 
        # est-ce que tu
        # peux aussi calculer la taille max d'un gène et modifier le 1000 du assert suivant stp ?

    def __getitem__(self, key) :
        
        return self.CONFIG[key]

    def __setitem__(self, key, value) :
        
        self.CONFIG[key] = value

    def __contains__(self, key) :
        
<<<<<<< HEAD
        # next_fitness = self.get_fitness(PARAMS["w_path_1"]+"params_seq.ini",PARAMS["w_path_1"])
        next_fitness = self.get_fitness(PARAMS["path_params_seq"]+"params_seq.ini")
=======
        return key in self.CONFIG

    def mutate(self, sim = 1, rep = 1) : 

        # if simulation is ended
        if self['SIMS_DONE'] >= self['N_SIM'] and \
           self['REPS_DONE'] >= self['N_REP'] :
                
                return

        sim = self.CONFIG['N_SIM'] if sim <= 0 else sim
        rep = self['N_REP'] if rep <= 0 else rep
    
        for repetition in range(self.rep, rep) :
            
            if repetition > 0 : self.restart_plasmid()
            
            for simulation in range(self.time, sim + 1) :
                
                print('%s [S: %d/%d R: %d/%d]'%(self['CONFIG_NAME'], 
                                                       self.time,
                                                       self['N_SIM'],
                                                       self.rep+1,
                                                       self['N_REP']))
            
                #MUTATION
                choice = np.random.choice(['Ins','Del','Inv'], p = self['PROBS']) 
                
                if LOG : print('\tOperation:\t%s'%(choice))
            
                apply_mut = self.do_my[choice]
                updated_data = apply_mut(self.data)

                #SELECTION
                next_fitness = self.get_fitness(updated_data)
                keep_new = self.keep_mutated(next_fitness)
                
                if keep_new :
                    
                    if LOG : print('\tPlasmid kept')
                    
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
                
                starts = (self.data['GFF']['seq'])[['start', 'end']].max()
                stops  = (self.data['GFF']['seq'])[['start', 'end']].min()
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
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
        
                self.save()
        
<<<<<<< HEAD
        if keep_new :
            print("\tmutate !")
            #UPDATE BEST
            # utils.save_data_to_path(updated_data,PARAMS["w_path_1"])
            self.data = copy.deepcopy(updated_data)
            self.fitness = next_fitness
            if choice == "Inv":
                self.genes = pd.concat([self.genes,self.genes_after_inversion])
            #UPDATE HISTORY
            else :
                self.append_plasmid_description()
        else :
            self.append_plasmid_description()

        self.hist_timestamp.append(self.time)
        self.hist_event.append(choice)
        self.hist_fitness.append(self.fitness)
        self.hist_kept.append(keep_new)
            ### RAJOUTER UNE DESCRIPTION COMPLETE DU GENOME.
        
    
    def get_fitness(self,params_file):
        proportions = simulation.start_transcribing_2(params_file,self.data) 
=======
    def get_fitness(self, data) :
        
        proportions = simulation.start_transcribing_2(self.config, data) 
        proportions = proportions/sum(proportions)
        
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        
        diff = -abs(proportions-self.target['expression'].values)/len(proportions)
        
        return(round(np.sum(diff),5))

    def keep_mutated(self,next_fitness):
        if LOG : print('\tFitness:\t%f -> %f'%(round(self.fitness,4), round(next_fitness,4)))
        if next_fitness>self.fitness:
            return(True)
        else : 
            alpha = math.exp(-self['ALPHA_C']*abs(self.fitness-next_fitness))
            if LOG : print('\tAlpha:    \t%f'%alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))
    
<<<<<<< HEAD
    def initialize_genes_description(self):
        start = self.data["GFF"]["seq"].ix[:,3].values
        end = self.data["GFF"]["seq"].ix[:,4].values
        orientation = self.data["GFF"]["seq"].ix[:,6].values
        length = abs(end-start)
        names = ["g"+str(i+1) for i in range(len(start))]

        e_start = np.copy(end)+1
        e_end = np.copy(np.roll(start,-1))-1
        e_end[-1] =  self.data["GFF"]["seq_length"]
        eee = list(zip(e_start,e_end))
        print("ici")
        print(eee)
        i = 0
        to_delete = []
        while i < len(eee):
            e_s, e_e = eee[i]
            for p in self.data["Prot"]["prot_pos"].values:
                if p>e_s and p<e_e:
                    eee.append((e_s,p-1))
                    eee.append((p+1,e_e))
                    to_delete.append(i)
                    break
            i+=1
        for i in sorted(to_delete,reverse=True):
            del eee[i]
        
        e_names = ["e"+str(i+1) for i in range(len(eee))]

        eee = np.array(eee)
        e_start = eee[:,0]

        e_end = eee[:,1]
        p = self.data["Prot"]["prot_pos"].values
        p_names = ["p"+str(i+1) for i in range(len(p))]

        names = names + e_names + p_names
        start = start.tolist() + e_start.tolist() + p.tolist()
        length = length.tolist() + abs(e_end-e_start).tolist() + [20 for i in range(len(p))] #changeable pour mise en valeur  
        orientation = orientation.tolist() + ["+" for i in range(len(eee))] + ["+" for i in range(len(p))]
        df = pd.DataFrame({"id" : names})
        ts = str(self.time)
        df["pos"+ts] = start
        df["length"+ts] = length
        df["dir"+ts] = orientation
        df.set_index("id",inplace=True)
        df = df.transpose()
        return(df)

    def append_plasmid_description(self):
        df = self.initialize_genes_description()
        self.genes = pd.concat([self.genes, df])
        return(0)
=======
    def normalize_plasmid(self):
    
        offset = self.data['TSS']['TSS_pos'][0] - 1

        self.data['GFF']['seq']['start'] -= offset
        self.data['GFF']['seq']['end'] -= offset
        self.data['Prot']['prot_pos'] = (self.data['Prot']['prot_pos']-offset)%self.data['GFF']['seq_length']

        self.data['TTS']['TTS_pos'] -= offset
        self.data['TSS']['TSS_pos'] -= offset

        return
    
    def check_config(self):

        min_gene_size = np.abs(np.min(self.data['TTS']['TTS_pos']-self.data['TSS']['TSS_pos']))
        
        assert self['U'] < min_gene_size
        assert np.all(self['PROBS'] >= 0.0)
        assert np.sum(self['PROBS']) == 1
        assert np.all(self.data['TSS']['TSS_pos'][0] > 0)
        assert np.all(self.data['TTS']['TTS_pos'][0] > 0)
        assert self.data['TSS']['TSS_pos'][0] == 1
        
        return
    
    def restart_plasmid(self) :
        
        self.rep += 1
        self.time = 0
        
        # Reading original data
        self.data = utils.read_plasmid_data(self.CONFIG)
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
        
        starts = (self.data['GFF']['seq'])[['start', 'end']].max()
        stops  = (self.data['GFF']['seq'])[['start', 'end']].min()
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

        if LOG : print('Plasmid configuration saved to\t%s'%location)

    def save_history(self):
    
        # Build history data frame
        
        subdict = { key : self.history[key] for key in self.history if key != 'plasmid' }
        
        df = pd.DataFrame(data = subdict)
    
        # saving location
        location = self['HISTORY_SAVE_PATH']
        
        # Save to CSV
        df.to_csv(path_or_buf = location, index = False, sep = '\t')
    
        if LOG : print('Fitness history saved to\t%s'%location)

    def save_config(self):
        
        self.config.set('PATHS', 'WPATH', str(self['WPATH']))
        self.config.set('STATE', 'SIMS_DONE', str(self.time-1))
        self.config.set('STATE', 'REPS_DONE', str(self.rep+1))
        
        file_name = self['WPATH'] + 'config.ini'
        
        with open(file_name, 'w') as config_file :
            
            self.config.write(config_file)
            
        if LOG : print('Config saved to\t%s'%file_name)
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f

    def U_inversion(self, data):
        
        updated_data = copy.deepcopy(data)
        
        l = data['GFF']['seq_length']     
        
        a = np.random.randint(0,l)
        b = np.random.randint(0,l)
        
        (a, b) = (b, a) if  a > b else  (a, b) # reversing if necessary

<<<<<<< HEAD
        start = updated_data["GFF"]["seq"].apply(lambda x : min(x["start"],x["end"]),axis=1)
        end = updated_data["GFF"]["seq"].apply(lambda x : max(x["start"],x["end"]),axis=1)
        # print(start,end)
=======
        start = updated_data['GFF']['seq'].apply(lambda x : min(x['start'],x['end']),axis=1)
        end = updated_data['GFF']['seq'].apply(lambda x : max(x['start'],x['end']),axis=1)
 
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
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
    
<<<<<<< HEAD
        updated_data['GFF']['seq'].columns = old_names
        data['GFF']['seq'].columns = old_names
        

        self.genes_after_inversion = self.genes.ix[-3:,:].copy()
        print(self.genes_after_inversion)
        real_positions = self.initialize_genes_description()
        #change e index only



=======
>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f
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
    
    def resume_plasmid(self):
        
        # Reading actual data
        self.data = utils.read_saved_plasmid(self.CONFIG)
        self.normalize_plasmid()
        
        self.rep = self['REPS_DONE'] - 1
        self.time = self['SIMS_DONE'] + 1
        
<<<<<<< HEAD
        path = PARAMS["w_path_1"]
        print("Saving evolution history to : "+path)
        
        self.history.to_csv(path_or_buf = path+"history.csv", index=False, sep=',',
                            columns=['time', 'event', 'kept', 'fitness'])
        

        self.genes = self.genes.transpose()
        print(self.genes)
        names = self.genes.index
        self.genes.reset_index(drop=True)
        self.genes.insert(0,"id",names.values)
        self.genes.to_csv(path_or_buf = path+"plasmid_description.csv", index=False, sep=',')
    
        # Save
    
        utils.save_data_to_path(self.data, PARAMS["w_path_1"] + 'final_')
=======
        # resume history
        df_hist = pd.read_table(self['WPATH'] + self['HISTORY_SAVE_NAME'])
        
        self.history['fitness'].extend(df_hist['fitness'])
        self.history['time'].extend(df_hist['time'])
        self.history['repetition'].extend(df_hist['repetition'])
        self.history['event'].extend(df_hist['event'])
        self.history['kept'].extend(df_hist['kept'])
        
        self.history['plasmid_size'].extend(df_hist['plasmid_size'])
        self.history['gene_ratio'].extend(df_hist['gene_ratio'])
        self.history['up_down_ratio'].extend(df_hist['up_down_ratio'])
        self.history['mean_space'].extend(df_hist['mean_space'])
        
        # resume fitness
        self.fitness = self.history['fitness'][-1]
        
        # resume plasmid
        plasmid_lines = ''
        
        with open( self['WPATH'] + self['PLASMID_SAVE_NAME']) as plasmid_file :
            
            plasmid_lines = plasmid_file.read().splitlines()
        
        self.history['plasmid'].extend(plasmid_lines[1:])
    
    def save(self) :
        
        self.save_history()
        self.save_plasmid_description()
        utils.save_data(self.data, self.CONFIG) #TODO save for each repetition !
        self.save_config()
        
        return

>>>>>>> f38c9629b6c9b1a56b829977c233834e70e8050f

        
        