

# A case for the existence of an evolution process driven by Indel and Inversion mutations

Note : This is a work in progress. 

We make use of the existing implementation of a transcription model taking into account
the physical structure of DNA and the coupling of neighboor genes, to proceed to the simulation
of the evolution process that could occur only trough the relative reorganisation of the genome.

Hypothesis : ... 

#### Authors

- Charles GAYDON
- Baptiste LAC

### Usage
This code is designed for Linux. 
First, use:

	conda env create --file EvoEnv.yml

to install a conda env - or install the dependencies manually using pip3.
Then activate the EvoEnv environment:

	source activate EvoEnv

 and run:

	python start_evol_simulation.py params_evo.ini 


#### Coding Conventions

indent = 4 spaces

Python version > 3.x