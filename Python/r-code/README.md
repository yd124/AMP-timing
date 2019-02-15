R scripts for processing python output
=================================


## sim2monolix.R

`sim2monolix.R` will convert sparse sampling python simulations (in wide format) into long-format for monolix. Will also merge on true times.

How to use (in bash) from within `AMP-timing/Python` working directory:

```
 Rscript r-code/sim2monolix.R sparse_data_filename.csv output_file_name.csv true_time_filename.csv
```

Working example from within `AMP-timing/Python` working directory (2019-02-14):

```
 Rscript r-code/sim2monolix.R sparsedata.csv sparsedata_mono.csv truetimes.csv
```

Only the sparsedata.csv is a required specification. `sim2monolix.R` will not overwrite output files!
