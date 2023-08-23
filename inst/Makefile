
dependencies:
	(cd .. && R --vanilla < inst/code/dependencies.R &)

intro:
	(cd .. && R --vanilla < inst/code/intro-binary.R &)

faithfulness:
	(cd .. && R --vanilla < inst/code/faithfulness.R &)

gcm-smooth:
	(cd .. && R --vanilla < inst/code/gcm-smooth.R &)

gcm-penalized:
	(cd .. && R --vanilla < inst/code/gcm-penalized.R &)

identifiability:
	(cd .. && R --vanilla < inst/code/identifiability.R &)

invariance-example:
	(cd .. && R --vanilla < inst/code/invariance-example.R &)

roc-invariance:
	(cd .. && R --vanilla < inst/code/roc-invariance.R &)

run-simulation:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 &)
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 2 &)
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 3 &)

vis-simulation:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 1 &)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 2 &)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 3 &)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 4 &)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 5 &)

all: intro faithfulness gcm-smooth gcm-penalized identifiability invariance-example

sim: roc-invariance run-simulation

vis: vis-simulation
