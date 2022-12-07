
# MicroMeta: Meta-analysis version of QCAT

I provide a meta-analysis framework for Microbiome data sets from
different studies based on the summary statistics generated by
Quasi-Conditional Association (QCAT) tests<sup>1</sup>. The
meta-analysis methods in this package include fixed effect tests, for
example, the FE-SKAT test and multivariate test, and random effect
tests, like, Het-SKAT and RE-SKAT test<sup>2</sup>. If the taxonomic
hierarchy is provided, the tests will be performed on the taxonomic tree
to localize the covariate-associated lineages.

*MicroMeta* comes with four datasets: data.meta, count.genus, meta and
tax. The data.meta dataset follows the structure of *QCAT_Meta* and
*QCAT_GEE_Meta* functions are designed. And data.meta\$OTU and
data.meta\$covariate are both lists contains 5 elements which are the
OTU counts and covariate of interest comes from different sources.
Meanwhile, data.meta, count.genus, meta and tax are the raw data sets of
data.meta.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# Plain installation
# install.packages("devtools")
devtools::install_github("DriWater/MicroMeta")

# For installation with vignette
devtools::install_github("DriWater/MicroMeta", build_vignettes = TRUE)
```

## Example

``` r
library(MicroMeta)
data("data.meta")
head(data.meta$Tax)
head(data.meta$OTU[[1]])
head(data.meta$covariate[[1]])
```

## QCAT_Meta

#### Fixed effect Meta-analysis method based on QCAT test with asymptotic p-value. (Tax Provided)

``` r
Method = "FE-MV"
one_part_MV <- QCAT_Meta(data.meta$OTU, data.meta$covariate, 1, Tax = data.meta$Tax, Method = Method, n.perm = NULL)
one_part_MV
```

#### Random effect Meta-analysis method based on QCAT test with permutation p-value performed on lineage.

``` r
Method = "Het-SKAT"
one_part_Het <- QCAT_Meta(data.meta$OTU, data.meta$covariate, 1, Tax = NULL, Method = Method, n.perm = 1000)
one_part_Het
```

For OTU table with excessive zero counts, we also employ the generalized
estimating equation (GEE) method to estimate the zero part of estimate
equations.

## QCAT_GEE_Meta

#### Fixed effect Meta-analysis method based on QCAT_GEE test with asymptotic p-value performed on lineage.

``` r
Method = "FE-VC"
zero_part_MV <- QCAT_GEE_Meta(data.meta$OTU, data.meta$covariate, 1, Tax = NULL, Method = Method, n.perm = NULL)
zero_part_MV 
```

#### Random effect Meta-analysis method based on QCAT_GEE test with permutation p-value. (Tax provided)

``` r
Method = "RE-SKAT"
zero_part_RE <- QCAT_GEE_Meta(data.meta$OTU, data.meta$covariate, 1, Tax = data.meta$Tax, Method = Method, n.perm = 200)
zero_part_RE
```

<div id="refs" class="references csl-bib-body">

<div id="ref-tang2017general" class="csl-entry">

<span class="csl-left-margin">1 </span><span
class="csl-right-inline">Tang Z-Z, Chen G, Alekseyenko AV, Li H. A
general framework for association analysis of microbial communities on a
taxonomic tree. *Bioinformatics* 2017; **33**: 1278–85.</span>

</div>

<div id="ref-MetaSKAT" class="csl-entry">

<span class="csl-left-margin">2 </span><span
class="csl-right-inline">Lee S. MetaSKAT: Meta analysis for SNP-set
(sequence) kernel association test. 2013.</span>

</div>

</div>
