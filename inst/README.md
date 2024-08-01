
# Replication materials

This folder contains the replication material and is structured as follows.

1. `code/` contains the `R` scripts for replicating all results in the
manuscript.

2. `helpers/` contains a helper package that is loaded in some of the 
replication scripts.

We provide a `Makefile` for conveniently running the scripts in `code/`
(code for how to manually reproduce the results can be found in the `Makefile`) 
as follows.

## Installing dependencies

Please run this first.

```
make dependencies
```

## Running the simulation studies

```
make sim-all
```

## Reproducing Figures in Section 4.4 and Appendix B2

```
make vis
```

## Reproducing all other results

```
make all
```

# Detailed list of files in `./code/`

* `case-study.R` -- SUPPORT2 data application
* `dependencies.R` -- Script for installing all dependencies
* `environment-parent-faithfulness-violated.R` -- Results in Appendix A5
* `faithfulness.R` -- Results in Appendix B5
* `gcm-penalized.R` and `gcm-smooth.R` -- Results in Appendix B6
* `identifiability.R` -- Results in Figure E1
* `intro-binary.R` -- Results in Figure 1
* `invariance-example.R` -- Results in Appendix D
* `nonparametric.R` -- Results in Appendix B7
* `run-simulation.R` -- Results in Appendix B (see `Makefile` for details)
* `vis-simulation.R` -- For visualizing the output of `run-simulation.R`
