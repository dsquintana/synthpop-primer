## Synthetic datasets: A non-technical primer for the biobehavioural sciences to promote reproducibility and hypothesis-generation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3520475.svg)](https://doi.org/10.5281/zenodo.3520475)

Synthetic datasets are an emerging method originally developed to permit the sharing of confidential census data. Synthetic datasets mimic real datasets by preserving their statistical properties and the relationships between variables. Importantly, this method also reduces disclosure risk to essentially nil as no record in the synthetic dataset represents a real individual. This is the accompanying R script for [my primer](https://psyarxiv.com/dmfb3/), which enables scholars to create synthetic datasets and assess their utility via the `synthpop` R package. By sharing synthetic datasets that mimic original datasets that could not otherwise be made open, researchers can ensure the reproducibility of their results and facilitate data exploration while maintaining participant privacy.

#### Run the analysis in your web browser

To launch a RStudio server instance and run my analysis scripts online, click [here](https://mybinder.org/v2/gh/dsquintana/synthpop-primer/master?urlpath=rstudio) or on on the "Launch Binder" badge below.

  <!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dsquintana/synthpop-primer/master?urlpath=rstudio)
  <!-- badges: end -->

Once the Rstudio server instance has loaded, run the commands in the "R_script.R" file.

Due to resource constraints of the RStudio server instance, the scripts that create Supplementary Figures 1-3 described in the [companion preprint](https://psyarxiv.com/dmfb3/) could not be included. These scripts can be found on the paper's [Open Science Framework page](https://osf.io/z524n/)

#### Run this analysis locally
To run the analysis locally in RStudio, download this repository as a zipped file. The R version and package versions are noted in the `sessionInfo.txt` file.
