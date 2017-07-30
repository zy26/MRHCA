# MRHCA

MRHCA is a non-parametric co-expression
analysis method for large association network: 

1.  it can be generally applied to association networks of most assessment methods; 
2.  it outputs exact significance level for each identified hub and module; 
3.  the method is with relatively small computational and memory consumption, hence can be applied to large association networks; 
4.  it is sensitive to the modules of weak associations; and  
5.  it enables overlapped modules in its outputs. 

Our analysis has demonstrated MRHCA can 
1.  deal with large association networks, 
2.  rigorously assess statistical significance for hubs and module sizes, 
3.  identify co-expression modules with low associations, 
4.  detect small and significant modules, and 
5.  allow genes to be present in more than one modules, compared with existing methods.

## Link to the related codes

The link to R and C++ codes which could reproduce the result of "[MRHCA: A nonparametric statistics based method for hub and co-expression module identification in large gene co-expression network](https://doi.org/10.1007/s40484-018-0131-z)".

-   [R codes](https://github.com/zy26/MRHCA-R)
-   [C++ codes](https://github.com/zy26/mrct)

## R Package Installation

To install the development version of MRHCA, you will need to install at
least the following packages from CRAN

``` r
install.packages("Rcpp")
```

For Windows users,
Rtools (https://cran.r-project.org/bin/windows/Rtools/) should also be
installed.

Then,

``` r
install.packages("devtools")
devtools::install_github("zy26/MRHCA")
```

## R Package Usage

``` r
file <- "https://github.com/zy26/mrct/raw/master/testdata/E_coli_anaerobic.txt"
x <- as.matrix(read.table(file, sep = ' ', header = TRUE, row.names = 1))
mr <- MRHCA::GetHubs(x)
```

## Fix the Results from C++ Version

If you have already obtained the results from the [C++ codes](https://github.com/zy26/mrct), you can further use the ```FixHubs``` function to filter and optimize the results. 

``` r
datafile <- "TCGA-COAD.htseq_fpkm.tsv"
file <- "TCGA-COAD.htseq_fpkm.tsv.id.txt"
emfile <- "TCGA-COAD.htseq_fpkm.tsv.txt"
mr <- MRHCA::FixHubs(datafile, file, emfile)
```