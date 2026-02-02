# ExSEnt: Extrema-Segmented Entropy Analysis of Time Series

Code accompanying:

ExSEnt: Extrema-Segmented Entropy Analysis of Time Series  
S. Kamali, F. Baroni, P. Varona  
DOI: https://doi.org/10.48550/arXiv.2509.07751

This repository implements ExSEnt, a set of entropy-based metrics computed on extrema-segmented time series. ExSEnt provides separate measures of (i) time irregularity, (ii) amplitude irregularity, and (iii) time–amplitude coupling irregularity, enabling a more granular characterization of dynamical structure than single-entropy summaries.

## Data sources and references

Daphnet dataset (Parkinson’s disease, freezing of gait; vertical ankle trajectory): https://doi.org/10.1109/TITB.2009.2036165  
EMG dataset (motor movement): https://doi.org/10.1093/gigascience/gix034  
Preprocessing pipeline description: https://doi.org/10.1016/j.compbiomed.2025.110874

Violin plot code (external dependency):  
Bechtold, Bastian (2016). Violin Plots for Matlab (GitHub).  
Repository: https://github.com/bastibe/Violinplot-Matlab  
DOI: https://doi.org/10.5281/zenodo.4559847

## Code overview (Matlab + Python)

1. compute_ExSEnt_metrics.m/compute_ExSEnt_metrics.py  
Main entry point for computing the three ExSEnt metrics: H_D, H_A, and H_DA (entropies of time, amplitude, and time–amplitude irregularities of the extrema-segmented time series). Corresponds to the Methodology section.  
Dependencies: extract_DA.m, sample_entropy.m and extract_DA.py, sample_entropy.py.

2. ExSEnt_for_Synthetic_signals.m  
Generates synthetic signals (Gaussian noise, pink noise, Brownian motion) and computes ExSEnt measures. Used for Fig. 2.

3. Logistic_ExEnt_bifur.m  
Computes and plots the logistic-map bifurcation diagram, color-coded by H_DA, and plots ExSEnt measures alongside SampEn.  
Related: Logistic_map_sample_r_plots.m generates sample time series and phase-space plots for three control-parameter values. Used for Fig. 3.

4. Rossler_Bifurcation_ExEnt.m  
Computes and plots the Rössler-system bifurcation diagram, color-coded by H_DA, and plots ExSEnt measures alongside SampEn.  
Related: Rossler_system_sample_Cvalue_plots.m generates sample time series and phase-space plots for four control-parameter values.  
Numerical integration: ode4.m. Used for Fig. 4.

5. Rulkov_ExEnt_bifur_sigma.m  
Computes and plots the Rulkov-map bifurcation diagram, color-coded by H_DA, and plots ExSEnt measures alongside SampEn.  
Related: Rulkov_sample_values_plot.m generates sample time series and phase-space plots for four control-parameter values.  
Related: Rulkov_Lagged_sample_phase_space_plots.m generates phase space and lagged phase space plots for four sigma values. Used for Figs. 5–6.

6. Daphnet dataset analysis (FOG vs NFG)  
This group of scripts extracts freezing-of-gait (FOG) and non-freezing (NFG) intervals and computes ExSEnt on band-passed trajectories, then generates violin plots comparing FOG vs NFG across two frequency bands.  
Daphnet_bandpassed_ExSEnt_violin_plots.m computes ExSEnt and generates violin plots.  
Daphnet_SamEn_violin_plots.m generates analogous violin plots for SampEn. Used for Figs. 9–10.

## Citation

If you use this code, please cite:

ExSEnt: Extrema-Segmented Entropy Analysis of Time Series  
S. Kamali, F. Baroni, P. Varona  
https://doi.org/10.48550/arXiv.2509.07751

## Rights

All rights reserved by the authors.
