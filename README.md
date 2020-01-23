# FAUST
### Full Annotation Using Shape-constrained Trees

The `FAUST` package implements the FAUST method described in [New interpretable machine learning method for single-cell data reveals correlates of clinical response to cancer immunotherapy](https://www.biorxiv.org/content/10.1101/702118v2).

The `FAUST` package requires `Rcpp` and `devtools`, and that a C++11 compiler is available.

## Installation

Currently `faust` must be installed from its source. It depends on the `scamp` package.

The most recent version can be installed from [github](https://github.com/FredHutch/faust) using [devtools](https://github.com/r-lib/devtools) in R. An installation (without vignettes) can be performed by first installing:

    tryCatch(installed.packages()["BiocManager","Version"],
	     error = function(e){
                install.packages("BiocManager")
             })	
    library(BiocManager)
    if (BiocManager::version()!='3.10') {
       BiocManager::install(version="3.10")
    }

    BiocManager::install("Biobase", update = FALSE, version = "3.10")
    BiocManager::install("flowCore", update = FALSE, version = "3.10")
    BiocManager::install("flowWorkspace", update = FALSE, version = "3.10")

    tryCatch(installed.packages()["devtools","Version"],
	     error = function(e){
                install.packages("devtools")
             })	
    library(devtools)
    devtools::install_github("RGLab/scamp")

Once these preliminary libraries are installed, run the command:

    devtools::install_github("RGLab/FAUST")
    
To build the vignettes during installation, instead run:

    tryCatch(installed.packages()["knitr","Version"],
	     error = function(e){
                install.packages("knitr")
             })	
    tryCatch(installed.packages()["rmarkdown","Version"],
	     error = function(e){
                install.packages("rmarkdown")
             })
    tryCatch(installed.packages()["ggdendro","Version"],
	     error = function(e){
                install.packages("ggdendro")
             })
    tryCatch(installed.packages()["remotes","Version"],
    	     error = function(e){
                install.packages("remotes")
             })	
    remotes::install_github("RGLab/FAUST", force = TRUE, build_vignettes = TRUE)

This takes longer since the vignettes must be built from source.

After loading `FAUST`, type `vignette('faustIntro')` to read a vignette discussing how to use the `FAUST` function in R.

## Citation

If you end up using `FAUST` to analyze cytometry datasets, please consider citing [New interpretable machine learning method for single-cell data reveals correlates of clinical response to cancer immunotherapy](https://www.biorxiv.org/content/10.1101/702118v2).