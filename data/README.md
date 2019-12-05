AMP Timing Data
===================

## RV217 datasets

All data set cleaning code is in `RV217-Analysis/data-cleaning`

- **RV217Master.csv**: raw RV217 data
- **RV217Clean.csv**: cleaned up raw RV217 data. Field names and classes standardized and cleaned up. Other munging. [Details here.](https://github.com/dbrvs/AMP-timing/blob/master/RV217-Analysis/data-cleaning/RV217-DataClean.md)
- **RV217Mono.csv**: Log VL (`log10VL`) by `day` (starting at 0 is first pos. VL) within first year of data (i.e., primary kinetics).
- **RV217MonoCells.csv**: Extension of **RV217Mono.csv**, contains multiple outcome variables for the cell subtypes.
- **RV217MonoCens.csv**: Extension of **RV217Mono.csv**, contains all negative VL observed before day 0 denoted as censored (`cens`).
- **RV217MonoVLPred.csv**: Extension of **RV217MonoCens.csv** including prediced viral loads (`log10VL_pred`) from APTIMA data (`APTIMA_num`). Two additional censor flags were added.  `cens_pred` does not consider predicted VLs as censored. `cens_firstpred` only considers the first predicted VL as uncensored (some IDs had multiple predicted VLs).