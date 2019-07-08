# FAUST
### Full Annotation Using Shape Constrained Trees

The `faust` package implements the FAUST data anlysis pipeline

## Installation

Currently `faust` must be installed from its source. It depends on the `scamp` (Selective Clustering Using Modes of Projection) clustering package. 

The most recent version can be installed from [github](https://github.com/FredHutch/faust) using [devtools](https://github.com/r-lib/devtools) in R. A quick installation (without vignettes) can be performed by:

    library(devtools)
    devtools::install_github("RGLab/scamp")
    devtools::install_github("RGLab/FAUST")
    
To build the vignettes during installation, run:

    library(devtools)
    devtools::install_github("RGLab/scamp")
    devtools::install_github("FredHutch/faust", build_vignettes=T)
    
This takes longer since the vignettes must be built from source.

After loading `faust`, type `vignette('faustIntro')` to read a vignette discussing how to use the `faust` function in R.
