ROW = 1

dependencies:
	(cd .. && Rscript --vanilla inst/code/dependencies.R > dependencies.out)

intro:
	(cd .. && Rscript --vanilla inst/code/intro-binary.R)

faithfulness:
	(cd .. && Rscript --vanilla inst/code/faithfulness.R)

gcm-smooth:
	(cd .. && Rscript --vanilla inst/code/gcm-smooth.R)

gcm-penalized:
	(cd .. && Rscript --vanilla inst/code/gcm-penalized.R)

nonparametric:
	(cd .. && Rscript --vanilla inst/code/nonparametric.R)

identifiability:
	(cd .. && Rscript --vanilla inst/code/identifiability.R)

invariance-example:
	(cd .. && Rscript --vanilla inst/code/invariance-example.R)

case-study:
	(cd .. && R CMD BATCH --vanilla inst/code/case-study.R)

sim-correct:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 10 10 100 0 $(ROW))

sim-link:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 2 10 10 100 0 $(ROW))

sim-hidden:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 3 10 10 100 0 $(ROW))

sim-roc:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 4 10 10 100 0 $(ROW))

sim-larger:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 5 10 10 100 0 $(ROW))

sim-test:
	(cd .. && Rscript --vanilla inst/code/run-simulation.R 1 2 2 2 1 $(ROW))

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

vis-larger:
	(cd .. && Rscript --vanilla inst/code/vis-simulation.R 6)

all: intro faithfulness gcm-smooth gcm-penalized nonparametric identifiability invariance-example case-study

vis: vis-main vis-app vis-hidden vis-link vis-larger

count:
	ls test.* | wc

correct:
	bash dispatch.sh 4000 run-correct.sh 10800

hidden:
	bash dispatch.sh 5600 run-hidden.sh 43200

link:
	bash dispatch.sh 4000 run-link.sh 10800

roc:
	bash dispatch.sh 500 run-roc.sh 7200

larger:
	bash dispatch.sh 5600 run-larger.sh 43200

test:
	bash dispatch.sh 2 run-test.sh 300

deps:
	bash dispatch.sh 1 run-install.sh 1200
