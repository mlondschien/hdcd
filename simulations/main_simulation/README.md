## README
This folder contains files required to replicate the main simulation results (Tables 1, 2) of [Change point detection for graphical models in presenceof missing values](https://arxiv.org/pdf/1907.05409.pdf). The file `simulation.R` needs to be called with 14 arguments `comment`, `n`, `p`, `Network`, `NA_method`, `optimizer`, `deletion`, `delete_fraction`, `seed`, `nrep`, `n_of_segments`, `segment_lengths`, e.g.
```
R simulation.R rfcd_9 500 100 RandomNetwork pair line mcar 0.1 2 100 4 70 120 120 190
```
The results of the simulation will be stored in an `rds` file.

The files `simulation_line.txt` and `simulation_section.txt` contain bash code to run simulations for different setups on different cores.

The `aggregate_simulation_results.R` file contains code to read in the `rds` files from the simulations and aggregate the simulation results.