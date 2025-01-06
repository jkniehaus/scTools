# seuratHelpR
#### A group of R functions to help manipulate and visualize scRNA-seq data (seurat objects) in R 

### description/origin
These functions are built for seurat versions >= 5. Many functions in this package are influenced by or originate entirely from Shekhar, K. et al. (https://doi.org/10.1016/j.cell.2016.07.054) and Loo L. et al. (https://doi.org/10.1038/s41467-018-08079-9)

### installation
```
library(devtools)

install_github('https://github.com/jkniehaus/seuratHelpR.git')
library(seuratHelpR)
```

Each R function is listed within the 'man' and information regarding utility and parameters can be found using `?function` in R.
e.g. `?binomcount.test`
