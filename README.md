# tree_scripts
Scripts to extract information from trees and tree distributions

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
