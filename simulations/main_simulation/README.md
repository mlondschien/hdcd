## README
This folder contains files required to replicate the main simulation results (Tables 1, 2) of [Change point detection for graphical models in the presence of missing values](https://arxiv.org/pdf/1907.05409.pdf). The file `simulation.R` needs to be called with 14 arguments `comment`, `n`, `p`, `Network` (one out of `"RandomNetwork"`, `"ChainNetwork"`), `NA_method` (one out of `"loh_wainwright_bias_correction"`, `"complete_observations"`, `"pairwise_covariance"` and `"average_imputation"`), `optimizer` (either `"line_search"` for BS or `"section_search"` for OBS), `deletion` (either `"mcar"` or `"blockwise"`), `delete_fraction`, `seed`, `nrep`, `n_of_segments`, `segment_lengths`, e.g.
```bash
R --vanilla --slave <simulation.R rfcd_9 500 100 RandomNetwork pairwise_imputation section_search mcar 0.1 2 100 4 70 120 120 190 >
```
Note the trailing whitespace.

The results of the simulation will be stored in an `rds` file, names via the parameters.

The files `simulation_line.txt` and `simulation_section.txt` contain bash code to run simulations for different setups.

The `aggregate_simulation_results.R` file contains code to read in the `rds` files from the simulations, aggregate the simulation results, and produce latex code for tables 1 and 2. Note that it might be necessary to adjust `directory`. The archive `/data/simulation_results.zip` contains results used to create the figures of the original document.
