# Measuring rank robustness in scored protein interaction networks

All work was carried out in R, using the following packages
* stringr
* igraph
* Matrix
* Brobdingnag
* parallel
* ggplot2
* ggthemes
* grid
* gridExtra

The workflow can be split in three parts: data preparation, metric extraction, and ranking and robustness analysis. Scripts for figure generation are also included.

## 1. Data preparation
Data preparation includes parsing interaction data, generating synthetic networks, and formatiing these into igraph objects and data frames so they can then be passed to metric extraction. The associated R script is
* `data_prep.R`

Note that while interaction scores are between 0 and 1, these are multiplied by 1000 and recorded and treated as integers. This is in keeping with how STRING data is stored.

### 1.1. Parsing data from STRING and HitPredict
Data files associated with this are:
* `string_yeast_data.RData`
* `string_ecoli_data.RData`
* `string_pvivax_data.RData`

### 1.2. Generating synthetic networks
Data files associated with this are:
* `pvx_syn_data.RData`
* `gnp_syn_data.RData`

## 2. Metric extraction
Node metric extraction uses the prepared data and results in data frames containing the measured raw node metrics for each network. Calculating LOUD natural connectivity is currently very intensive, and if you're planning on using the code on large networks, you might need to remove it. The associated R scripts are:
* `metric_extraction_main.R`
* `metric_extraction_aux.R`

The auxiliary file contains standard and LOUD node metric extraction functions, and the main contains an example of how these are applied to the SYN-PVX network. An example of the split raw output is given in 
* `pvx_syn_A.RData`
* `pvx_syn_B.RData`
* `pvx_syn_C.RData`

These are generated in `metric_extraction_main.R`, and are then merged to a data frame `pvx.syn.metrics.df`, which is stored in 
`F_pvx_syn_metrics.RData`. Dataframes for the other networks follow the same naming convention of `NETWORK.NAME.metrics.df` and are stored in their respective `F_NETWORKNAME.RData` files.

## 3. Ranking and robustness analysis
All ranking and robustness analysis is handled by
* `robustness_analysis_main.R`
* `ronbustness_analysis_aux.R`

Like with metric extraction, the main file contains a worked example for the SYN-PVX network, and the auxiliary file contains all relevant functions. Output of all steps are stored in the `F_NETWORKNAME.RData` files, along with the raw metric values.

## 4. Figure generation
Plots are made using `ggplot`. The theme, colour scheme, labellers, etc. are in `plotPreamble.R`. Figure and table generation itself is in `figure_generation.R`. The relevant data across all six networks for rank continuity, identifiability, and instability can be found in
* `continuityPlots.RData`
* `identifiabilityPlots.RData`
* `instabilityPlots.RData`
