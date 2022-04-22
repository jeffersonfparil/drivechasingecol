# Slow and steady wins the race: spatial and stochastic processes and the failure of suppression gene drives

Code, and data repository.

## Authors:
- Jeff Paril ([ORCID: 0000-0002-5693-4123](https://orcid.org/0000-0002-5693-4123))
- Ben Phillips ([ORCID: 0000-0003-2580-2336](https://orcid.org/0000-0003-2580-2336))

## Quickstart guide:
- Edit `src/slurmer.sh` to match the configuration of your computing cluster.
- Run `src/driveChaseEcolRun.R` via `src/slurmer.sh`.
- Once all simulation runs have finished, edit `src/driveChaseEcolAnalysis.R` to match the set-up of your machine.
- Run `driveChaseEcolAnalysis.R` and assess the results.

## Code organisation:
- `driveChaseEcolRun.R` contains the main execution function for running each simulation.
- `driveChaseEcolFunctions.R` contains the growth, reproduction, dispersal, drive conversion, and other related fuctions called by the main simulation function.
- `PointMetricsDriveChaseEcol2D.c` contains the density and spatial coordinates calculation functions.
- `slurmer.sh` is a shell script used to generate multiple shell scripts to run on a multiple compute nodes on a computing cluster.
- `driveChaseEcolAnalysis.R` is the statistical analysis and plotting script.

## Notes:
- Output figures of the analysis are in `out/no_TADS-NstaR5/`.
- `out/Figure_4-recursive_function_dataset.rds` is generated using the ouput of `driveChaseEcolFunctions.R::evolve()` function at $\sigma=10$, and the recursive equations are explicitly implemented in `src/misc/simple_plot_timeseries_ribbon.R`, which was used to generate `out/no_TADS-Nstar5/Figure_4-timeseries-simulated-deterministic.svg`.
