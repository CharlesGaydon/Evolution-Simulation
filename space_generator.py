import os
import configparser
import pandas as pd
import numpy as np


def read_base(path):
    
    config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
    # to preserve capital letters
    config.optionxform = str
    config.read(path)
    
    print('%s read as base for generator'%path)
    
    return config
    
def save_probs(config, probs, prefix = 'PART', output_dir='paramsfiles/'):
    
    name = prefix + '_%d'
    path = name + '.ini'
    files = []
    
    if not os.path.isdir(os.path.join(os.getcwd(), output_dir)) :
        os.system('mkdir ' + output_dir)
    
    for i, prob in enumerate(probs) :
        
        a,b,c = prob
        
        config.set('SIMULATION', 'PINV', str(a))
        config.set('SIMULATION', 'PDEL', str(b))
        config.set('SIMULATION', 'PINS', str(c))
        
        config.set('INPUTS', 'CONFIG_NAME', name%i)
    
        output = output_dir + path%i
    
        with open(output, 'w') as config_file :
            
            config.write(config_file)
    
        files.append(output)
    
        print('Generated %s with probs %s'%(output, np.round(np.array(prob), 2)))
    
    return files
    
def save_alphas(config, alphas, prefix = 'APART', output_dir='paramsfiles/'):
    
    name = prefix + '_%d'
    path = name + '.ini'
    files = []
    
    if not os.path.isdir(os.path.join(os.getcwd(), output_dir)) :
        os.system('mkdir ' + output_dir)
    
    for i, alpha in enumerate(alphas) :
        
        
        config.set('SIMULATION', 'ALPHA_C', str(int(alpha)))
        
        config.set('INPUTS', 'CONFIG_NAME', name%i)
    
        output = output_dir + path%i
    
        with open(output, 'w') as config_file :
            
            config.write(config_file)
    
        files.append(output)
    
        print('Generated %s with alpha %s'%(output, alpha))
    
    return files
    
def simplex_generator(nsteps):
    
    slist = []
    
    x = np.linspace(0,1,nsteps)
    
    a,b,c = np.meshgrid(x,x,x)
    
    for i in np.linspace(0,1,nsteps) :
        for j in np.linspace(0,1,nsteps) :
                
                v = i, j, 1-(i+j)
                
                if np.all(np.array(v) >= 0) :
                
                    slist.append(v)
    
    print('Generated %d points'%len(slist))
    
    return slist
 
def save_as_task(files, file_name, output_dir = ''):
    
    if not os.path.isdir(os.path.join(os.getcwd(), output_dir)) :
        os.system('mkdir ' + output_dir)
    
    with open(output_dir + file_name, 'w') as mfile:
        
        mfile.write('\n'.join(files))
        
    print('Saved task as %s'%(output_dir+file_name))
    
    return
    
def main(args):
    
    path = args[1]
    comp = int(args[2])
    
    config = read_base(path)
    points = simplex_generator(comp)
    files = save_probs(config, points, output_dir='paramsfiles/light_%d/'%comp,
    prefix='LIGHT')
    save_as_task(files, 'light_%d.txt'%comp)
    
    #config = read_base(path)
    #points = np.linspace(200,1000,comp)
    #files = save_alphas(config, points, output_dir='paramsfiles/aparts_%d/'%comp)
    #save_as_task(files, 'ALPH_%d.txt'%comp)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

