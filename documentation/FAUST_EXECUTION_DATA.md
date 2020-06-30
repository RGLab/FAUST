![faust_logo](images/logos/faust_logo.png)

# Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

-   [Execution Data](#execution-data)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Execution Data

This is where you can get a high level understanding of the data that `FAUST` creats and uses during run-time

# Directory Structure - AKA `faustData`

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
        └── SAMPLE 03
            └── ...
```

## `exhaustiveFaustCountMatrix.rds`

TODO

## `expUnitData`

TODO

### `EXPERIMENTAL UNIT N`

TODO

## `faustCountMatrix.rds`

TODO

## `gateData`

TODO

## `metaData`

TODO

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

## `sampleData`

TODO
