import os

scripts = []

scripts += ['graphs.r']
scripts += ['rep.r']

for s in scripts :
    
    os.system('Rscript %s'%s)
    