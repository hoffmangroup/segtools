#!/usr/bin/env Rscript
list_packages <- c("cluster", "plyr", "grDevices", "latticeExtra", "reshape2", "lattice", "grid", "Cairo", "RColorBrewer")
missing_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) install.packages(missing_packages,repos='http://cran.us.r-project.org')
