# Install all dependencies for reproducing the results in the paper
# LK 2023

### Install remotes for installing all other packages
install.packages("remotes", repos = "https://cloud.r-project.org")

### GitHub packages
remotes::install_github("igraph/rigraph")
remotes::install_github("strobl/RCIT")

### Bioconductor packages
remotes::install_bioc(c("graph", "Rgraphviz", "RBGL"))

### Install all other packages from CRAN
pkgs <- c("MASS", "SimDesign", "coin", "colorspace", "cotram", "generalhoslem",
          "glmnet", "multcomp", "pROC", "patchwork", "tidyverse", "tram", "mlt",
          "ranger", "beepr", "pcalg")
remotes::install_cran(pkgs)


### Install tramicp and helpers locally
remotes::install_local(force = TRUE)
remotes::install_local("inst/helpers", force = TRUE)
