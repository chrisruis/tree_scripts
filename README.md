# tree_scripts
Scripts to extract information from trees and tree distributions

## bootstrap_TempEst_rttd_date.R
Calculates the significance of a collection date vs root-to-tip correlation using bootstrapping

Before running this script, import the tree into TempEst, set the best fitting root if needed and use "Export Data" to save a text file containing sample collection dates and root-to-tip distances. This script runs on this text file

The correlation between sample collection date and root-to-tip distance is calculated with the real sample collection dates. The significance of this correlation is calculated by randomly assigning the dates to tips and recalculating the correlation a given number of times. The proportion of these randomisations that result in a correlation at least as strong as with the real data is reported as the p-value

This script takes the text file from TempEst and the number of bootstraps to carried out. 1000 bootstraps is recommended. The script prints the R squared correlation with the real dates and the p-value on this correlation

To run:

RScript bootstrap_TempEst_rttd_date.R TempEst_text_file.txt number_of_bootstraps

E.g.

RScript bootstrap_TempEst_rttd_date.R rttd_date.txt 1000

## population_increase_distribution_BEAST.py
Calculates a distribution of increase in relative genetic diversity from a posterior distribution

This uses information within the Bayesian skyline population model parameters and will not work for other population models

Takes the .trees and .log files from BEAST. Identifies the first increase in relative genetic diversity from the PopSizes columns in the log file and identifies the date of this increase using the corresponding number of nodes in the GroupSizes columns and the node heights in the respective tree

Prints the proportion of trees in the posterior distribution that support an increase in relative genetic diversity. Writes the date of the first increase in relative genetic diversity for each MCMC step in which an increase is inferred

The first increase in relative genetic diversity is identified as the first window in the Bayesian skyline plot whose PopSize is more than x% above the baseline PopSize. The value of x% is set with option -p and is 100 by default. A value of 100 means the PopSize needs to double above baseline to be inferred as an increase

To convert node heights in the tree to dates, the date of the latest sequence in the tree needs to be supplied in decimal format (e.g. 2015.54) with -d

By default, the script expects the log file to be in BEAST2 format, in which case the PopSize and GroupSize column names should contain PopSize and GroupSize, respectively. This will not be the case with BEAST1 output. If using BEAST1 files, use option -b 1 which will switch so the script expects the PopSize and GroupSize columns to contain popSize and groupSize, respectively

To run:

python3 population_increase_distribution_BEAST.py -t BEAST.trees -l BEAST.log -d latest_sample_date -o output_file_name.txt

## population_change_support_BEAST.py
Calculates the proportion of sampled MCMC steps that support a change in relative genetic diversity within a given time window

This uses information within the Bayesian skyline population model parameters and will not work for other population models

Takes the .trees and .log files from BEAST and a window of interest. Use the dates in the tree and the GroupSizes in the log file to identify dates at which the relative genetic diversity changes. For each tree, checks if there is a change in relative genetic diversity within the window based on these dates and the PopSizes. This change is defined relative to the PopSize at the start of the window of interest

Prints the proportion of trees in the posterior distribution that support a change in relative genetic diversity within the window of interest

Supply the window of interest with -w. This takes 2 decimal numbers separated by a space, e.g. using "-w 1990 2000" will test for a change between 1990 and 2000

By default, looks for an increase in relative genetic diversity. To instead look for a decrease use --decrease

To convert node heights in the tree to dates, the date of the latest sequence in the tree needs to be supplied in decimal format (e.g. 2015.54) with -d

By default, the script expects the log file to be in BEAST2 format, in which case the PopSize and GroupSize column names should contain PopSize and GroupSize, respectively. This will not be the case with BEAST1 output. If using BEAST1 files, use option -b 1 which will switch so the script expects the PopSize and GroupSize columns to contain popSize and groupSize, respectively

To run:

python3 population_change_support_BEAST.py -t BEAST.trees -l BEAST.log -d latest_sample_date -w window_start window_end
