# FWI Synchronicity

This repository contains all the codes and data used for the Fire Weather Index (FWI) with observations (past and present) and perturbed (future). 
All codes are available and organized for reproducibility.

## Data

The data from the perturbed FWI simulations used in El Garroussi (2024) are available in a public repository at Zenodo. You can access them through the following link:

[https://doi.org/10.5281/zenodo.10458186](https://doi.org/10.5281/zenodo.10458186)

The observed FWI data are available at:

[https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.0e89c522?tab=overview](https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.0e89c522?tab=overview)

## Reference
El Garroussi, S. (2024). 30-Year Canadian Fire Weather Index Simulations over Europe: CMIP6-Informed Temperature and Precipitation Perturbations [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10458186


## Scripts
This repository contains all the scripts used for the synchronicity analysis in our study:

### For the "past" and "present" periods (Table 1, Figure 2 and Figures S1-S4):
#### 1. A1_synchronicity_reg_season_p95
Performs regional synchronicity analysis for the season, considering the 95th percentile.
#### 2. A2_synchronicity_reg_season_50_38
Carries out regional synchronicity analysis for the season, considering FWI thresholds 50 and 38.
#### 3. A3_Tab1_synchronicity
Generates a table showing significant changes in all analyzed parameters. Corresponds to Table 1.

### For the "past" and "future" periods (Figure 3 and Figures S6-S10):
#### 4. B1_synchronicity_perturbed_50_38
Performs a perturbation analysis to study the robustness of synchronicity results considering FWI thresholds 50 and 38.
#### 5. B2_synchronicity_reg_season_p95
Executes a specific perturbation analysis for the 95th percentile.
#### 6. C1_plot_trend_perturbed
Generates trend plots to visualize synchronicity in different perturbation contexts.
#### 7. C2_plot_surface_perturbed
Generates additional surface plots for synchronicity analysis under perturbed conditions.

### To obtain the seasonal 95th percentile (Figure S1):
#### 8. D1_p95_season

### To analyze trends in the Fire Weather Index (Figure S5):
#### 9. E1_trend_FWI

## Script Usage

Each script is designed to be executed in the R environment. Ensure that the necessary packages are installed before running the scripts. You can run each script as follows:

```R
# Example of how to run a script in R
source("path/to/script/A1_synchronicity_reg_season_p95.R")
