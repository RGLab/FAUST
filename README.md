# FAUST
### Full Annotation Using Shape-constrained Trees

The `FAUST` package implements the FAUST method described in a forthcoming manuscript as an R package.

It assumes `Rcpp` is installed, and that a C++11 compiler is available.

Currently, `FAUST` works only on MAC OSX and Linux.  

## Installation

Currently `faust` must be installed from its source. It depends on the `scamp` package.

The most recent version can be installed from [github](https://github.com/FredHutch/faust) using [devtools](https://github.com/r-lib/devtools) in R. A quick installation (without vignettes) can be performed by:

    library(devtools)
    devtools::install_github("RGLab/scamp")
    devtools::install_github("RGLab/FAUST")
    
To build the vignettes during installation, run:

    library(devtools)
    devtools::install_github("RGLab/scamp")
    devtools::install_github("RGLab/FAUST", build_vignettes=T)
    
This takes longer since the vignettes must be built from source.

After loading `FAUST`, type `vignette('faustIntro')` to read a vignette discussing how to use the `FAUST` function in R.
