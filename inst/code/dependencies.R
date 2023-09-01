# Install all dependencies for reproducing the results in the paper
# LK 2023

.libPaths(c("~/tutu/lib", .libPaths()))

### Install remotes for installing all other packages
install.packages("remotes", repos = "https://cloud.r-project.org")

### Install all other packages from CRAN
pkgs <- c("MASS", "SimDesign", "coin", "colorspace", "cotram", "generalhoslem",
          "glmnet", "multcomp", "pROC", "patchwork", "tidyverse", "tram", "mlt",
          "ranger", "beepr", "pcalg")
remotes::install_cran(pkgs)

### Install tramicp and helpers locally
remotes::install_local()
remotes::install_local("inst/helpers")
