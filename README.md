# FAUST

### Full Annotation Using Shape-constrained Trees

The `FAUST` package implements the FAUST method described in [New interpretable machine learning method for single-cell data reveals correlates of clinical response to cancer immunotherapy](https://www.biorxiv.org/content/10.1101/702118v2).

Full Annotation Using Shape-constrained Trees (FAUST) is a method to perform discovery and annotation of phenotypes in single-cell data from flow and mass cytometry experiments.

This implementation targets the gating set data structure in the [`flowWorkspace`](https://bioconductor.org/packages/release/bioc/html/flowWorkspace.html) Bioconductor package.

When applied to a cytometry study stored in a gating set, the `faust` R function creates a directory called "faustData" located at the setting of a parameter called "projectPath".

The principal output is an **annotated count matrix** `faustCountMatrix.rds` which is written to the "faustData" directory.

- **Rows** in the count matrix correspond to **samples** collected in the single cell experiment.
- The **columns** in the count matrix are **cell populations** discovered by the pipeline. 
- The **columns** are *annotated* by a selected subset of *markers* used in conducting the experiment.
- The **column annotations** define, in terms of these *markers*, the **phenotypes** of all **cell populations** discovered by the pipeline.
- **Entries** in the count matrix correspond the **number of cells** in a sample that belong to a discovered cell population. 

This count matrix can be loaded into R using the `readRDS` function.

## Installation

The `FAUST` package requires `Rcpp` and `devtools`, and that a C++11 compiler is available.

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