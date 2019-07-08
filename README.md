# FAUST
### Full Annotation Using Shape-constrained Trees

The `FAUST` package implements the FAUST method described in a forthcoming manuscript as an R package.

The `FAUST` package requires `Rcpp` and `devtools`, and that a C++11 compiler is available.

Currently, `FAUST` works only on MAC OSX and Linux.  

## Installation

Currently `faust` must be installed from its source. It depends on the `scamp` package.

The most recent version can be installed from [github](https://github.com/FredHutch/faust) using [devtools](https://github.com/r-lib/devtools) in R. A quick installation (without vignettes) can be performed by:

    tryCatch(installed.packages()["BiocManager","Version"],
	     error = function(e){
                install.packages("BiocManager")
             })	
    library(BiocManager)
    BiocManager::install("Biobase", update = FALSE)
    BiocManager::install("flowCore", update = FALSE)
    BiocManager::install("flowWorkspace", update = FALSE)
    BiocManager::install("flowWorkspaceData", update = FALSE)

    library(devtools)
    devtools::install_github("RGLab/scamp")
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
    devtools::install_github("RGLab/FAUST", build_vignettes=T)
    
This takes longer since the vignettes must be built from source.

After loading `FAUST`, type `vignette('faustIntro')` to read a vignette discussing how to use the `FAUST` function in R.
