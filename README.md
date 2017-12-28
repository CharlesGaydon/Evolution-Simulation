

# A case for the existence of an evolution process driven by Indel and Inversion mutations

We make use of the existing implementation of a transcription model taking into account
the physical structure of DNA and the coupling of neighboor genes, to proceed to the simulation
of the evolution process that could occur only trough the relative reorganisation of the genome.

Hypothesis : ...

#### Authors

- Charles GAYDON
- Baptiste LAC

### Organisation

- `display/` contains the R-script used to display the fitness graph;
- `docs/` contains the server used to display the plasmid evolution;
- `environments/` contains the target environment's description for the expression of the genes;
- `paramsfiles/` contains the `.ini` files used to describe the simulation's parameters;
- `plasmids/` contains the different plasmids used in the simulations;
- `simulations/` contains the simulation's results (history and plasmids).


### Usage

This code is designed for Linux. 

1. First, install a conda env containing the dependencies

	`conda env create --file EvoEnv.yml`


2. Then, activate the EvoEnv environment

	`source activate EvoEnv`

3. Run a simulation
	
	`python start_evol_simulation.py paramsfiles/params_evo.ini`

### Fitness graphs

1. Install dependencies in R

	`install.packages(ggplot2)`

2. Set-up the correct working directory

	`setwd("your_path")`

3. Copy a `history.csv` file next to `display/graphs.r`

4. Run the script 

	`Rscript graphs.r`

5. Enjoy the result

![a fitness graph](https://github.com/CharlesGaydon/Evolution-Simulation/blob/master/display/example_1.pdf)

#### Coding Conventions

- indent = 4 spaces
- Python version > 3.x
