# Install all dependencies for reproducing the results in the paper
# LK 2023

### Install remotes for installing all other packages
install.packages("remotes", repos = "https://cloud.r-project.org")

### GitHub packages
remotes::install_github("igraph/rigraph")
remotes::install_github("ericstrobl/RCIT")

### Bioconductor packages
remotes::install_bioc(c("graph", "Rgraphviz", "RBGL"))

### Install all other packages from CRAN
pkgs <- c("MASS", "SimDesign", "coin", "colorspace", "cotram", "generalhoslem",
          "glmnet", "multcomp", "pROC", "patchwork", "tidyverse", "tram", "mlt",
          "ranger", "beepr", "pcalg", "tramME")
remotes::install_cran(pkgs)


### Install tramicp and helpers locally
remotes::install_local(force = TRUE)
remotes::install_local("inst/helpers", force = TRUE)

### Create directories for outputs
if (!dir.exists(tp <- file.path("inst", "figures")))
  dir.create(tp, recursive = TRUE)

if (!dir.exists(tp <- file.path("inst", "results")))
  dir.create(tp, recursive = TRUE)
