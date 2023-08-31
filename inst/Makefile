ROW = 1

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
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 20 20 300 0 $(ROW))

sim-link:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 2 20 20 300 0 $(ROW))

sim-hidden:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 3 20 20 300 0 $(ROW))

sim-roc:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 4 20 20 300 0 $(ROW))

sim-test:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 2 2 2 1 $(ROW))

sim-gcm:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 20 4 100 1 $(ROW))

sim-all: sim-correct sim-link sim-hidden sim-roc

vis-main:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 1)

vis-app:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 2)

vis-hidden:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 3)

vis-link:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 4)

vis-wald:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 5)

all: intro faithfulness gcm-smooth gcm-penalized identifiability invariance-example

vis: vis-main vis-app vis-hidden vis-link vis-wald

count:
	ls test.* | wc

correct:
	bash dispatch.sh 12000 run-correct.sh

hidden:
	bash dispatch.sh 12000 run-hidden.sh

link:
	bash dispatch.sh 12000 run-link.sh

roc:
	bash dispatch.sh 1500 run-roc.sh

test:
	bash dispatch.sh 2 run-test.sh
