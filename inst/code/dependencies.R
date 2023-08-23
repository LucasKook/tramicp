# Install all dependencies for reproducing the results in the paper
# LK 2023

### Install remotes for installing all other packages
install.packages("remotes", repos = "https://cloud.r-project.org")

### Install tramicp and helpers locally
remotes::install_local(dependencies = TRUE)
remotes::install_local("inst/helpers", dependencies = TRUE)

### Install RCIT from GitHub
remotes::install_github("ericstrobl/RCIT")

### Install all other packages from CRAN
pkgs <- c("MASS", "SimDesign", "coin", "colorspace", "cotram", "future.apply",
          "generalhoslem", "glmnet", "multcomp", "pROC", "patchwork", "tramME",
          "tramnet", "tidyverse")
remotes::install_cran(pkgs)
