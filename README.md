PAAQD2: A code-refined version of PAAQD
===============================================

Installation
------------------------

-   Install the latest version from [GitHub](https://github.com/masato-ogishi/PAAQD2) as follows:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/PAAQD2")
```

-   You might be prompted to install some packages before installling Repitope. Follow the message(s).

Loading
------------------

``` r
library(PAAQD2)
```

Usage
-----------------------------------
``` r
PAAQD2(peptideSet=c("YILDLQPEN","MCLRFLSKI","GLSPAITKY","WWFQSSMSK","ATATELNNA","VVVVVVVVV"))
```

Reference
-----------------------------------

Saethang, T., Hirose, O., Kimkong, I., Tran, V. A., Dang, X. T., Nguyen, L. A. T., ... & Satou, K. (2013). PAAQD: Predicting immunogenicity of MHC class I binding peptides using amino acid pairwise contact potentials and quantum topological molecular similarity descriptors. Journal of immunological methods, 387(1), 293-302.
