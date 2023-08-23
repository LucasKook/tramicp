
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

sim-correct:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1)

sim-link:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 2)

sim-hidden:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 3)

sim-roc:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 4)

sim-all: sim-correct sim-link sim-hidden sim-roc

vis-simulation:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 1)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 2)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 3)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 4)
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 5)

all: intro faithfulness gcm-smooth gcm-penalized identifiability invariance-example

vis: vis-simulation
