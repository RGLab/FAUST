![faust_logo](images/logos/faust_logo.png)

# Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

-   [FAUST Inputs](#faust-inputs)
    -   [Method Signature](#method-signature)
-   [Required Inputs](#required-inputs)
    -   [gatingSet](#gatingset)
        -   [Example](#example)
            -   [Using a CytoBank Gating Set](#using-a-cytobank-gating-set)
            -   [Using Raw FCS Files](#using-raw-fcs-files)
            -   [Using a FlowJo Gating Set](#using-a-flowjo-gating-set)
            -   [Using a Legacy FlowWorkspace Gating Set](#using-a-legacy-flowworkspace-gating-set)
            -   [Using a FlowWorkspace Gating Set](#using-a-flowworkspace-gating-set)
    -   [startingCellPop](#startingcellpop)
        -   [activeChannels](#activechannels)
            -   [Manual Creation](#manual-creation)
            -   [Loading From a File](#loading-from-a-file)
    -   [channelBounds](#channelbounds)
    -   [imputationHierarchy](#imputationhierarchy)
    -   [experimentalUnit](#experimentalunit)
    -   [projectPath](#projectpath)
    -   [depthScoreThreshold](#depthscorethreshold)
    -   [selectionQuantile](#selectionquantile)
    -   [nameOccuranceNum](#nameoccurancenum)
    -   [supervisedList](#supervisedlist)
    -   [debugFlag](#debugflag)
    -   [threadNum](#threadnum)
    -   [seedValue](#seedvalue)
    -   [drawAnnotationHistograms](#drawannotationhistograms)
    -   [annotationsApproved](#annotationsapproved)
    -   [densitySubSampleThreshold](#densitysubsamplethreshold)
    -   [densitySubSampleSize](#densitysubsamplesize)
    -   [densitySubSampleIterations](#densitysubsampleiterations)
    -   [archDescriptionList](#archdescriptionlist)
    -   [plottingDevice Optional Inputs](#plottingdevice-optional-inputs)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# FAUST Inputs

This document explains this implementation of FAUST's input requirements.

For this document we are using the [faust(...) method's](../R/faust.R#L217) signature to determine the inputs.

## Method Signature

```R
faust <- function(gatingSet,
                  startingCellPop,
                  activeChannels=flowWorkspace::markernames(gatingSet),
                  channelBounds="",
                  experimentalUnit="",
                  imputationHierarchy="",
                  projectPath=normalizePath("."),
                  depthScoreThreshold=0.01,
                  selectionQuantile=0.50,
                  nameOccuranceNum=0,
                  supervisedList=NA,
                  debugFlag=FALSE,
                  threadNum=1,
                  seedValue=123,
                  drawAnnotationHistograms=TRUE,
                  annotationsApproved=FALSE,
                  densitySubSampleThreshold=1e6,
                  densitySubSampleSize=1e6,
                  densitySubSampleIterations=1,
                  archDescriptionList=
                      list(
                          targetArch=c("singleCPU")
                      ),
                  plottingDevice="pdf"
                  )
{
    ...
}
```

# Required Inputs

## gatingSet

This is a [FlowWorkspace Gating set](https://github.com/RGLab/flowWorkspace). FAUST will work with legacy and modern gating sets from that library.

### Example

#### Using a CytoBank Gating Set

Please see the [Cytoverse Example Page](https://cytoverse.org/examples/index.html)

#### Using Raw FCS Files

TODO

#### Using a FlowJo Gating Set

Please see the [Cytoverse Example Page](https://cytoverse.org/examples/index.html)

#### Using a Legacy FlowWorkspace Gating Set

```R
library(flowWorkspace)

path_to_gating_set = file.path("/", "Users", "example_username", "Desktop", "example_gating_set")
path_to_converted_gating_set = file.path("/", "Users", "example_username", "Desktop", "converted_example_gating_set")

flowWorkspace::convert_legacy_gs(from=path_to_gating_set, to=path_to_converted_gating_set)

gating_set_to_use <- flowWorkspace::load_gs(path_to_converted_gating_set)

faust <- function(gating_set_to_use, ...)
```

#### Using a FlowWorkspace Gating Set

```R
library(flowWorkspace)

path_to_gating_set = file.path("/", "Users", "example_username", "Desktop", "example_gating_set")

gating_set_to_use <- flowWorkspace::load_gs(file_path_to_gating_set)

faust <- function(gating_set_to_use, ...)
```

## startingCellPop

This is the `gate` that `FAUST` should start at. If you're unsure set this to be `root` gate

### activeChannels

This is a file that specifies the active channels of the gating set to use when running `FAUST`.

⚠️**Warning**⚠️ - This MUST be the same length as the `channelBounds` that are provided as well.

#### Manual Creation

```R
active_channels_to_use <- c("CD14", "CD3E", "CD4", "IL2RA", "NCAM1", "CCR7", "CD19", "PTPRC", "CD8A")

faust <- function(active_channels_to_use, ...)
```

#### Loading From a File

```R
active_channels_to_use <- c("CD14", "CD3E", "CD4", "IL2RA", "NCAM1", "CCR7", "CD19", "PTPRC", "CD8A")
active_channels_file_path <- file.path("/", "Users", "example_username", "Desktop", "example_active_channels.rds")

saveRDS(active_channels_to_use, file= active_channels_file_path)

loaded_active_channels_file_path <- loadRDS(active_channels_file_path)

faust <- function(loaded_active_channels_file_path, ...)
```

## channelBounds

⚠️**Warning**⚠️ - This MUST be the same length as the `activeChannels` that are provided as well.

## imputationHierarchy

## experimentalUnit

## projectPath

## depthScoreThreshold

## selectionQuantile

## nameOccuranceNum

## supervisedList

## debugFlag

## threadNum

## seedValue

## drawAnnotationHistograms

## annotationsApproved

## densitySubSampleThreshold

## densitySubSampleSize

## densitySubSampleIterations

## archDescriptionList

## plottingDevice Optional Inputs
