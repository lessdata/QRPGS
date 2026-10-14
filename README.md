# Quantile regression-based polygenic score (QR-PGS)

This repository provides code and documentation for polygenic prediction and predictive uncertainty quantification using quantile regression, as described in our paper, "Uncertainty-aware polygenic prediction reveals individual variation in phenotype predictability and distributional risk profiles".



### Install dependent packages in R

- data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)
- dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)
- tidyr [https://CRAN.R-project.org/package=tidyr](https://CRAN.R-project.org/package=tidyr)
- stringr [https://CRAN.R-project.org/package=stringr](https://CRAN.R-project.org/package=stringr)
<!-- - splines2 [https://CRAN.R-project.org/package=splines2](https://CRAN.R-project.org/package=splines2) -->
<!-- - ggplot2 (for downstream analyses) [https://CRAN.R-project.org/package=ggplot2](https://CRAN.R-project.org/package=ggplot2) -->
- fsQRPPA (for model training) [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/)



### Model training with user's own phenotype and genotype data

Users can train their own models using the [fsQRPPA](https://arxiv.org/abs/2601.02826) method. The fsQRPPA R package, available at [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/), is required for model training.

An example, [train_qrpgs.R](train_qrpgs.R), is provided to demonstrate the model develpment.



## Citation

TBA.


