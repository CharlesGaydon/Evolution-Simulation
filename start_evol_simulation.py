# coding: utf-8
import sys
import evol_simulation as evol

if len(sys.argv)==1:
    print("Usage : python start_evol_simulation.py params_file.ini")
    exit(0)
elif len(sys.argv)==2:
    INI_file = sys.argv[1]
    evol.start_evol_simulation(INI_file)
else : 
    print("Usage : python start_evol_simulation.py params_file.ini")
    exit(0)

