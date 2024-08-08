# Organization
This file introduces the sub-directories `raw_data`, `code`, `estimation` and `figures` for producing the empirical study results in our submission. The figures in our main paper and supplement for empirical study are reproduced and saved in directory `empirical_study/figures` and our reported p-values in main paper are computed in R scripts under `code` directory. We will specify these details further when introducing the `code` section.

## raw data

The directory `empirical_study/raw_data/` contains our raw data obtained from the website https://www.metoffice.gov.uk/research/climate/maps-and-data/historic-station-data; Based on data integrity, we selected data from 27 stations. The website provides the names and location information of our selected stations, which we have copied into the `locations.csv` file and Table B.5 of our supplement.


## code

We recommend running codes under directory `empirical_study/code/` in following order:

### Results in main paper
1. Running `deseason.R` to pre-process the raw data by deseasonalization and missing-value interpolation as described in main paper. The R script would output `data_processed.csv` under directory `empirical_study/estimation/`.
2. Running `estimation_trend.R` to generate monotone estimation results upon the preprocessed data. The R script would output `high-dim_2023.RData` under directory `empirical_study/estimation/`.
3. Running `SCB_high_feasible.R` based on the estimation results to apply the bootstrap Algorithm 1 in our main paper. The R script would output `bootstrap_samples_2023.csv` under directory `empirical_study/estimation/`.
4. Running `test_results.R` to reproduce the Figure 2 in the Section 8 of main paper. Moreover, the p-values reported in the Section 8 would be printed by the R script.
5. Running `Figure-1.R` to reproduce the Figure 1 in main paper; Running `monotonicity.R` to print out the p-values of the Chetverikov's and Bowman's monotonicity test as described in page 2 of the main paper. 
6. For the results in supplement, one can run `reg_estimaton.R` to generate estimations for regression model in the Section B.4 of supplement. The R script would output `regression_sun.RData` under directory `empirical_study/estimation/`. Then running `reg_SCB.R` to obtain bootstrap sample `boot_reg.RData` in `empirical_study/estimation/`. Finally running `Figure-B.3.R` to reproduce the Figure B.3 in supplement as well as the p-values.























