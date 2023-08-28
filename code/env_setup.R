library(renv)

renv::init(bioconductor = TRUE)

renv::install("tidyverse")
renv::install("data.tree")
renv::install("furrr")
renv::install("Matrix")
renv::install("viridis")
renv::install("seriation")
renv::install("igraph")
renv::install('vroom')
renv::install("bioc::GenomicRanges")
renv::install("bioc::AnnotationHub")

renv::snapshot()
