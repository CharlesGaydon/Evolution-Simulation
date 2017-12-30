

# A case for the existence of an evolution process driven by Indel and Inversion mutations

We make use of the existing implementation of a transcription model taking into account
the physical structure of DNA and the coupling of neighboor genes, to proceed to the simulation
of the evolution process that could occur only trough the relative reorganisation of the genome.

# To-Do List

[x] implement simulation recover/resume
[x] implement simulation repetition
[x] implement multithreading
[x] correct statistic (mean_space)
[x] implement regular saving
[ ] correct plasmid saving for each repetition (**not relevant I think**)
[ ] implement automated graph generation
[ ] find other statistics
[ ] handle rare exception in simulation computing... It wastes time

# Hypothesis

- unit is 150 wide (insertions/deletions)
- the smallest gene is larger than the unit
- ...

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

3. Start a new simulation
	
	`python start_evol_simulation.py paramsfiles/standard.ini`
	
4. Or resume a simulation

	`python start_evol_simulation.py simulations/STD_XXXXXXXX/config.ini`

### Fitness graphs

1. Install dependencies in R

	`install.packages('ggplot2')`
	`install.packages('ggtern')`
	`install.packages('plotly')`

2. Set-up the correct working directory _if needed_

	`setwd("your_path")`

3. Copy a `history.csv` file next to `display/graphs.r`

4. Run the script 

	`Rscript graphs.r`

5. Enjoy the results

![a fitness graph](https://github.com/CharlesGaydon/Evolution-Simulation/blob/master/display/example_graph.png)

#### Coding Conventions

- indent = 4 spaces
- Python version > 3.x
