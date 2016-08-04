## Analysis of Data

### Data Files
OTU picking was done with Ninja at 97% identity against the Green Genes database.
Alpha diversity was calculated with alpha_diversity.py in QIIME v 1.8.0

### Dependencies

The analysis is done using R. The following R packages are required:

* ape
* vegan
* beeswarm
* RColorBrewer
* reshape2
* plyr
* grid
* gridExtra
* ggplot2

### Running the code
You can download the repository and enter the main directory. Using R, you can call the run.r script, which will run the full analysis. Do so from the command line with:

```
Rscript bin/run.r

```
or using R studio set your directory to the main repository folder (where the README.md is) and then type:

```
source("bin/run.r")
```
 