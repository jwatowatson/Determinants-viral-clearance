# Determinants-viral-clearance
Statistical analysis of an individual patient data meta-analysis for temporal changes in viral clearance kinetics and optimal design of the PLATCOV trial. The PLATCOV trial is registered at clinicaltrials.gov number [NCT05041907](https://clinicaltrials.gov/study/NCT05041907).

This work is licensed under a Creative Commons Attribution 4.0 International License.

[![alt text](https://upload.wikimedia.org/wikipedia/commons/e/e1/CC_BY_icon.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview
This github repo provides the data and code for the statistical analysis of the article "Temporal changes in SARS-CoV-2 clearance kinetics and the optimal design of antiviral pharmacodynamic studies: an individual patient data meta-analysis of a randomised, controlled, adaptive platform study (PLATCOV)" published in The Lancet Infectious Diseases. 

The analysis in this repo included:
1. The meta-analysis for individual patient data in unblinded patients of the PLATCOV trials between September, 2021 and October, 2023 to describe temporal changes in SARS-CoV-2 clearance kinetics. Temporal trends of baseline viral densities and viral clearance rates are described using a [penalized B-spline fit](https://github.com/milkha/Splines_in_Stan).
   * The analysis is described in _01_temporal_splines_analysis.qmd_ and _01_run_local_temporal_splines_analysis.R_.
   * Data used in this analysis is _Unblinded_meta_analysis.csv_.

2. The boostrapping analysis for the impacts follow-up duration on the z scores.
   * The analysis is described in _02_bootstrapping_analysis.qmd_ and _02_run_local_bootstraps_analysis.R_.
   * Data used in this analysis are:
     * _Ineffective_analysis.csv_
     * _Paxlovid_Molnupiravir_analysis.csv_
     * _Paxlovid_recent_analysis.csv_
     * _REGN_analysis.csv_
     * _Remdesivir_analysis.csv_
3. The boostrapping analysis for the impacts sampling frequencies on the z scores.
   * The analysis is described in _03_bootstrapping_analysis_sampling_fq.qmd_ and _03_run_local_bootstraps_analysis_sampling_fq.R_.
   * Data used in this analysis are:
     * _Ineffective_analysis.csv_
     * _Paxlovid_Molnupiravir_analysis.csv_
     * _Paxlovid_recent_analysis.csv_
     * _REGN_analysis.csv_
     * _Remdesivir_analysis.csv_

Stan models used in this analysis are in the _Stan_models_ folder:
* _Temporal_spline_additive.stan_
* _Temporal_spline_proportional.stan_: This is the main model for this analysis.

Any questions or comments or if any bugs spotted drop me a message at **jwatowatson@gmail.com** or **phrutsamon.wongnak@gmail.com**. 
