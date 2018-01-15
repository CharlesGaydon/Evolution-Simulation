

# A case for the existence of an evolution process driven by Indel and Inversion mutations

We make use of the existing implementation of a transcription model taking into account
the physical structure of DNA and the coupling of neighboor genes, to proceed to the simulation
of the evolution process that could occur only trough the relative reorganisation of the genome.

# To-Do List

- [x] implement simulation recover/resume
- [x] implement simulation repetition
- [x] implement multithreading
- [x] correct statistic (mean_space)
- [x] implement regular saving
- [x] correct plasmid saving for each repetition (**not relevant I think**)
- [x] implement automated graph generation
- [x] find other statistics
- [x] handle rare exception in simulation computing... It wastes time

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
- `simulations/` contains the simulation's results (history and plasmids);
- `tasks/` contains tasks to launch pooled simulations.

### Usage

This code is designed for Linux. 

1. First, install dependencies

	`conda env create --file EvoEnv.yml`

2. Then, activate the development environment

	`source activate EvoEnv`

3. Start a new simulation
	
	`python start_evol_simulation.py paramsfiles/standard.ini`
	
4. Or start new simulation**s** from a task file on 4 processsors
	
	`python start_evol_simulation.py @tasks/tasks.txt --nproc 4`
	
4. Or resume a stopped simulation

	`python start_evol_simulation.py simulations/XXX_XXXXXXXX/config.ini`

5. Use `python start_evol_simulation.py --help` for further information

6. Graphs are automatically generated in the `simulations/` folder !
