
# Replication materials

This folder contains the replication material and is structured as follows.

1. `helpers/` contains a helper package that is loaded in some of the 
replication scripts.

2. `code/` contains the `R` scripts for replicating all results in the
manuscript.

3. `results/` contains the saved results from the simulation study (which
takes several days to run on a 20-core machine), in order to reproduce the
figures without having to re-run the simulations. (This folder has not been
checked into GitHub, but a link will be made available shortly.)

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
