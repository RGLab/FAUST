![faust_logo](images/logos/faust_logo.png)

# Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

-   [Overview](#overview)
-   [`faustData`](#faustdata)
    -   [`exhaustiveFaustCountMatrix.rds`](#exhaustivefaustcountmatrixrds)
    -   [`expUnitData`](#expunitdata)
        -   [`EXPERIMENTAL UNIT N`](#experimental-unit-n)
            -   [`expUnitExprs.rds`](#expunitexprsrds)
            -   [`expUnitRes.rds`](#expunitresrds)
            -   [`expUnitToSampleLookup.rds`](#expunittosamplelookuprds)
            -   [`scampClusterLabels.rds`](#scampclusterlabelsrds)
            -   [`scampExpUnitComplete.rds`](#scampexpunitcompleterds)
    -   [`faustCountMatrix.rds`](#faustcountmatrixrds)
    -   [`gateData`](#gatedata)
        -   [`root_rawGateList.rds`](#root_rawgatelistrds)
        -   [`root_resList.rds`](#root_reslistrds)
        -   [`root_resListPrep.rds`](#root_reslistpreprds)
        -   [`root_selectedChannels.rds`](#root_selectedchannelsrds)
    -   [`metaData`](#metadata)
    -   [`plotData`](#plotdata)
        -   [`hist_MARKER_ab_1.pdf`](#hist_marker_ab_1pdf)
        -   [`histograms`](#histograms)
            -   [`MARKER #`](#marker-)
            -   [`scoreLines.pdf`](#scorelinespdf)
    -   [`sampleData`](#sampledata)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Overview

This is where you can get a high level understanding of the data that `FAUST` creats and uses during run-time. This data is stored in the `faustData` directory and if needed can be used for a more rigorous understanding of what has occurred during processing.

# `faustData`

In its current form, `FAUST` will create a working directory that contains all of the data used for clustering.

This directory structure looks like this.

```bash
└── faustData
    ├── exhaustiveFaustCountMatrix.rds
    ├── expUnitData
    │   ├── EXPERIMENTAL UNIT 01
    │   │   └── ...
    │   ├── EXPERIMENTAL UNIT 02
    │   │   ├── expUnitExprs.rds
    │   │   ├── expUnitRes.rds
    │   │   ├── expUnitToSampleLookup.rds
    │   │   ├── scampClusterLabels.rds
    │   │   └── scampExpUnitComplete.rds
    │   └── EXPERIMENTAL UNIT N
    │       └── ...
    ├── faustCountMatrix.rds
    ├── gateData
    │   ├── root_rawGateList.rds
    │   ├── root_resList.rds
    │   ├── root_resListPrep.rds
    │   └── root_selectedChannels.rds
    ├── metaData
    │   ├── activeChannels.rds
    │   ├── allScampClusterNames.rds
    │   ├── analysisMap.rds
    │   ├── channelBounds.rds
    │   ├── channelBoundsUsedByFAUST.rds
    │   ├── colNameMap.rds
    │   ├── depthMat.rds
    │   ├── firstALReady.rds
    │   ├── forceList.rds
    │   ├── initSelC.rds
    │   ├── madeResMats.rds
    │   ├── parsedGS.rds
    │   ├── phenotypeElbowValue.rds
    │   ├── possibilityList.rds
    │   ├── preferenceList.rds
    │   ├── sanitizedCellPopStr.rds
    │   ├── scampClusterNames.rds
    │   ├── scampNameSummary.rds
    │   ├── scoreMat.rds
    │   ├── selectionList.rds
    │   └── startingCellPop.rds
    ├── plotData
    │   ├── hist_MARKER_ab_1.pdf
    │   ├── histograms
    │   │   ├── MARKER 01
    │   │   │   └── ...
    │   │   ├── MARKER 02
    │   │   │   ├── EXPERIMENTAL UNIT 01.pdf
    │   │   │   ├── ...
    │   │   │   └── EXPERIMENTAL UNIT N.pdf
    │   │   └── MARKER N
    │   │       └── ...
    │   ├── scampNamesPlot.pdf
    │   └── scoreLines.pdf
    └── sampleData
        ├── SAMPLE 01
        │   └── ...
        ├── SAMPLE 02
        │   ├── annotationMatrix.csv
        │   ├── exhaustiveFaustAnnotation.csv
        │   ├── exprsMat.rds
        │   ├── faustAnnotation.csv
        │   ├── resMat.rds
        │   └── scampAnnotation.csv
        └── SAMPLE N
            └── ...
```

## `exhaustiveFaustCountMatrix.rds`

TODO

---

## `expUnitData`

TODO

### `EXPERIMENTAL UNIT N`

TODO

#### `expUnitExprs.rds`

TODO

#### `expUnitRes.rds`

TODO

#### `expUnitToSampleLookup.rds`

TODO

#### `scampClusterLabels.rds`

TODO

#### `scampExpUnitComplete.rds`

TODO

---

## `faustCountMatrix.rds`

TODO

---

## `gateData`

TODO

### `root_rawGateList.rds`

TODO

### `root_resList.rds`

TODO

### `root_resListPrep.rds`

TODO

### `root_selectedChannels.rds`

TODO

---

## `metaData`

TODO

---

## `plotData`

TODO

### `hist_MARKER_ab_1.pdf`

TODO

### `histograms`

TODO

#### `MARKER #`

TODO

#### `scoreLines.pdf`

TODO

---

## `sampleData`

TODO
