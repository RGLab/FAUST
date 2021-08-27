# FAUST

### Full Annotation Using Shape-constrained Trees

The `FAUST` package implements the FAUST method described in [New interpretable machine learning method for single-cell data reveals correlates of clinical response to cancer immunotherapy](https://www.biorxiv.org/content/10.1101/702118v2).

The `FAUST` package requires `Rcpp` and `devtools`, and that a C++11 compiler is available.

## Installation

Currently `faust` must be installed from its source. It depends on the `scamp` package. 

The most recent version can be installed from [github](https://github.com/FredHutch/faust) using [devtools](https://github.com/r-lib/devtools) in R. An installation (without vignettes) can be performed by first installing:

```
    tryCatch(installed.packages()["BiocManager","Version"],
	     error = function(e){
                install.packages("BiocManager")
             })	
    library(BiocManager)
    if (BiocManager::version()!='3.13') {
       BiocManager::install(version="3.13")
    }

    BiocManager::install("Biobase", update = FALSE)
    BiocManager::install("flowCore", update = FALSE)
    BiocManager::install("flowWorkspace", update = FALSE)

    tryCatch(installed.packages()["devtools","Version"],
	     error = function(e){
                install.packages("devtools")
             })	
    library(devtools)
    devtools::install_github("RGLab/scamp")
```

Once these preliminary libraries are installed, run the command:

```
    devtools::install_github("RGLab/FAUST")
```

To build the vignettes during installation, instead run:

```
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
```

This takes longer since the vignettes must be built from source.

After loading `FAUST`, type `vignette('faustIntro')` to read a vignette discussing how to use the `FAUST` function in R.

## MacOS Installation Note 

MacOS users should make sure developer tools and XCode are installed. If you run into installation issues, the function `pkgbuild::check_build_tools()` may be useful in debugging them. Additionally, the information about [R compiler tools for MacOS found here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) may also be helpful.

## Linux Installation Note 

If you are trying to install FAUST on Linux with GCC version greater than 8, and you encounter installation errors, the following note may help resolve these errors: [protocol buffer issue link](https://github.com/protocolbuffers/protobuf/issues/5144#issuecomment-688723405).

## Citation

If you end up using `FAUST` to analyze cytometry datasets, please consider citing [New interpretable machine learning method for single-cell data reveals correlates of clinical response to cancer immunotherapy](https://www.biorxiv.org/content/10.1101/702118v2).