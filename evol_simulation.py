# coding: utf-8
import pandas as pd
import numpy as np
import copy
import os
import simulation
import utils
import math
import time

LOG = False

'''
Authors : Charles GAYDON, Baptiste LAC
Oct 2017- February 2018
'''

def start_evol_simulation(config_path) :

        # create plasmid with config
        config = utils.read_config(config_path)

        p1 = Plasmid(config)

        try:
                p1.mutate(0,0)
                os.system('python display/run_scripts.py %s'%('' + p1['WPATH']))
                return(-1, p1['WPATH'])
        except:
                print('------ AN ERROR OCCURED -----')
                print('\tPlasmid name: %s'%p1['CONFIG_NAME'])
                os.system('python display/run_scripts.py %s'%('' + p1['WPATH']))
                return(1, p1['WPATH'])
                
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
        self.flag_resume = False
        self.do_my = {'Ins':self.U_insertion,
                      'Del':self.U_deletion,
                      'Inv':self.U_inversion } #TODO : add inversion

        self.config = config
        self.CONFIG = utils.parse_config(self.config)

        assert np.sum(self['PROBS']) != 0

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
            
            self.flag_resume = True
        
            self.resume_plasmid()

        self.check_config()

        print('Plasmid     | %s'%self['PLASMID_PATH'])
        print('Environment | %s'%self['TARGET_PATH'])
        print('Parameters  | %s'%self['CONFIG_NAME'])

    def __getitem__(self, key) :

        return self.CONFIG[key]

    def __setitem__(self, key, value) :

        self.CONFIG[key] = value

    def __contains__(self, key) :

        return key in self.CONFIG

    def mutate(self, sim = 1, rep = 1) :

        # if simulation is ended
        if self['SIMS_DONE'] >= self['N_SIM'] and \
           self['REPS_DONE'] >= self['N_REP'] :
                return

        sim = self.CONFIG['N_SIM'] if sim <= 0 else sim
        rep = self['N_REP'] if rep <= 0 else rep

        for repetition in range(self.rep, rep) :

            if repetition > 0 and not self.flag_resume:
                
                 self.restart_plasmid()
                 self.flag_resume = False

            for simulation in range(self.time, sim + 1) :

                #MUTATION
                choice = np.random.choice(['Ins','Del','Inv'], p = self['PROBS'])

                print('%s [S: %d/%d R: %d/%d]\t%s'%(self['CONFIG_NAME'],
                                       self.time,
                                       self['N_SIM'],
                                       self.rep + 1,
                                       self['N_REP'],
                                       choice))


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

                gene_ratio = ((self.data['TTS']['TTS_pos'] - \
                                    self.data['TSS']['TSS_pos']).abs()).sum()/plasmid_size

                UD_ratio = (self.data['TSS']['TUorient'] == '+').sum() / \
                           (self.data['TSS']['TUorient'] == '-').sum() if \
                           (self.data['TSS']['TUorient'] == '-').sum() != 0 else \
                           -1

                starts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
                stops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)
                starts_V = (stops[:-1] + 1).reset_index(drop=True)
                stops_V  = (starts[1:] - 1).reset_index(drop=True)
                lengths_V = stops_V - starts_V + 1
                lengths_V = np.concatenate(lengths_V , self.data['GFF']['seq_length'] - stops.iat[-1])
                mean_space = lengths_V.mean()

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

                self.update_plasmid_description()

                self.save()
                
                self.time += 1
                
                pass
            
            pass
        
        return

    def get_fitness(self, data) :

        proportions = simulation.start_transcribing_2(self.config, data)
   
        proportions = proportions/proportions.sum()

        #FONCTIONNE EN L'ABSENCE D'INVERSION DE GENES -> adapter le calcul qd les genes changent d'ordre ans gff
        return -(proportions-self.target['expression']).abs().sum()/len(proportions)

    def keep_mutated(self,next_fitness):

        if LOG : print('\tFitness:\t%f -> %f'%(round(self.fitness,4), round(next_fitness,4)))

        if next_fitness > self.fitness:
            return(True)
        else :
            alpha = math.exp(-self['ALPHA_C']*np.abs(self.fitness-next_fitness))
            if LOG : print('\tAlpha:    \t%f'%alpha)
            return(np.random.choice([True,False], p = [alpha,1-alpha]))

    # WILL BE DELETED
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
    # WILL BE DELETED
    def append_plasmid_description(self):
        df = self.initialize_genes_description()
        self.genes = pd.concat([self.genes, df])
        return(0)

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

        gene_ratio = ((self.data['TTS']['TTS_pos'] - self.data['TSS']['TSS_pos']).abs()).sum() \
                        /plasmid_size

        UD_ratio = (self.data['TSS']['TUorient'] == '+').sum() / \
                   (self.data['TSS']['TUorient'] == '-').sum() if \
                   (self.data['TSS']['TUorient'] == '-').sum() != 0 else \
                   -1

        starts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
        stops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)
        starts_V = (stops[:-1] + 1).reset_index(drop=True)
        stops_V  = (starts[1:] - 1).reset_index(drop=True)
        lengths_V = stops_V - starts_V + 1
        lengths_V = np.concatenate(lengths_V , self.data['GFF']['seq_length'] - stops.values[-1])
        mean_space = lengths_V.mean()

        self.history['fitness'].append(self.fitness)
        self.history['time'].append(0)
        self.history['repetition'].append(self.rep)
        self.history['event'].append('Beg')
        self.history['kept'].append(True)

        self.history['plasmid_size'].append(plasmid_size)
        self.history['gene_ratio'].append(gene_ratio)
        self.history['up_down_ratio'].append(UD_ratio)
        self.history['mean_space'].append(mean_space)

        self.update_plasmid_description()

        self.time += 1

        return

    def update_plasmid_description(self):

        tim = self.time
        typ = 'G'
        rep = self.rep

        starts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
        stops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)
        lengths = stops - starts + 1

        for index in range(len(starts)) :

            pos = starts.iat[index]
            ori = self.data['TSS']['TUorient'].iat[index]
            lth = lengths.iat[index]

            ori = (1 if ori == '+' else -1)

            self.history['plasmid'] += ['%d\t%d\t%s\t%d\t%d\t%d'%(rep, tim, typ, pos, lth, ori)]


        typ = 'V'

        starts_V = (stops.iloc[:-1] + 1).reset_index(drop=True)
        stops_V =  (starts.iloc[1:] - 1).reset_index(drop=True)
        lengths_V = stops_V - starts_V + 1

        ori = 0

        for index in range(len(starts_V)) :

            pos = starts_V.iat[index]
            lth = lengths_V.iat[index]

            self.history['plasmid'] += ['%d\t%d\t%s\t%d\t%d\t%d'%(rep, tim, typ, pos, lth, ori)]


        self.history['plasmid'] += ['%d\t%d\t%s\t%d\t%d\t%d'%(rep, tim, typ,
                                                              stops.iat[-1] + 1,
                self.data['GFF']['seq_length'] - stops.iat[-1], ori)]

        typ = 'P'
        ori = 0
        lth = 1

        for index in range(len(self.data['Prot'])) :

            pos = self.data['Prot']['prot_pos'].iat[index]

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

    def U_inversion(self, data):

        updated_data = copy.deepcopy(data)

        l = data['GFF']['seq_length']

        a = np.random.randint(1,l)
        b = np.random.randint(1,l)

        (a, b) = (b, a) if  a > b else  (a, b) # reversing if necessary

        rstarts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
        rstops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)

        while   (np.logical_or(
                                np.logical_and(a >= rstarts, a<=rstops),
                                np.logical_and(b >= rstarts, b<=rstops))).any() \
                or \
                (np.logical_or(a == data['Prot']['prot_pos'],
                               b == data['Prot']['prot_pos'])).any() :

                #print('rejected cuts %d / %d'%(a,b))

                a = np.random.randint(1,l)
                b = np.random.randint(1,l)

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
        updated_data['GFF']['seq'].sort_values(by='start', inplace=True)
        updated_data['TSS'].sort_values(by='TSS_pos', inplace=True)
        updated_data['TTS'].sort_values(by='TTS_pos', inplace=True)
        updated_data['Prot'].sort_values(by='prot_pos', inplace=True)

        return(updated_data)

    def U_deletion(self, data) :

        updated_data = copy.deepcopy(data)

        l = data['GFF']['seq_length']

        #localisation
        s1 = np.random.randint(1, l + 1 - self['U'] + 1) #entre [1 et 29851] par exemple
        s2 = s1 + self['U'] - 1                          #par exemple start = 29851, stop = 30000 (donc 150 bases)

        rstarts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
        rstops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)

        #sinon relancer
        while np.logical_and(s1 >= rstarts, s1 <= rstops).any() or \
              np.logical_and(s2 >= rstarts, s2 <= rstops).any() or \
              np.logical_and(s1 >= self.data['Prot']['prot_pos'],  \
                             s2 <= self.data['Prot']['prot_pos']).any() :

            s1 = np.random.randint(1, l + 1 - self['U'] + 1)
            s2  = s2 + self['U'] - 1


        #deletion
        updated_data['TTS']['TTS_pos']-= (self.data['TTS']['TTS_pos']>s2)*self['U']
        updated_data['TSS']['TSS_pos']-= (self.data['TSS']['TSS_pos']>s2)*self['U']
        updated_data['Prot']['prot_pos'] -= (self.data['Prot']['prot_pos']>s2)*self['U']
        updated_data['GFF']['seq']['start'] -= (self.data['GFF']['seq']['start']>s2)*self['U']
        updated_data['GFF']['seq']['end'] -= (self.data['GFF']['seq']['start']>s2)*self['U']

        updated_data['GFF']['seq_length'] = l - self['U']

        return(updated_data)

    def U_insertion(self,data) :

        updated_data = copy.deepcopy(data)

        l = data['GFF']['seq_length']

        #localisation
        s = np.random.randint(1, l+1)

        rstarts = (self.data['GFF']['seq'])[['start', 'end']].min(axis=1).reset_index(drop=True)
        rstops  = (self.data['GFF']['seq'])[['start', 'end']].max(axis=1).reset_index(drop=True)


        while np.logical_and(s >= rstarts, s <= rstops).any() or \
              (s == data['Prot']['prot_pos']).any() :

            s = np.random.randint(1, l+1)

        #insertion
        updated_data['TTS']['TTS_pos']+= (self.data['TTS']['TTS_pos']>s)*self['U']
        updated_data['TSS']['TSS_pos']+= (self.data['TSS']['TSS_pos']>s)*self['U']
        updated_data['Prot']['prot_pos'] += (self.data['Prot']['prot_pos']>s)*self['U']

        updated_data['GFF']['seq']['start'] += (self.data['GFF']['seq']['start']>s)*self['U']
        updated_data['GFF']['seq']['end'] += (self.data['GFF']['seq']['start']>s)*self['U']

        updated_data['GFF']['seq_length'] += self['U']

        return(updated_data)

    def resume_plasmid(self):

        # Reading actual data
        self.data = utils.read_saved_plasmid(self.CONFIG)
        self.normalize_plasmid()

        self.rep = self['REPS_DONE'] - 1
        self.time = self['SIMS_DONE']

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

        # WILL BE DELETED
        # <<<<<<< HEAD
        #         path = PARAMS["w_path_1"]

        #         self.genes = self.genes.transpose()
        #         names = self.genes.index
        #         self.genes.reset_index(drop=True)
        #         self.genes.insert(0,"id",names.values)
        #         self.genes.to_csv(path_or_buf = path+"plasmid_description.csv", index=False, sep=',')

        return




