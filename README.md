
# MicroMeta: Meta-analysis version of QCAT

I provide a meta-analysis framework for Microbiome data sets from
different studies based on the summary statistics generated by
Quasi-Conditional Association (QCAT) tests<sup>1</sup>. The
meta-analysis methods in this package include fixed effect tests, for
example, the FE-SKAT test and multivariate test, and random effect
tests, like, Het-SKAT and RE-SKAT test<sup>2</sup>. If the taxonomic
hierarchy is provided, the tests will be performed on the taxonomic tree
to localize the covariate-associated lineages.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DriWater/MicroMeta")
```

## To do list:

I plan to finish the code for the zero-part model (QCAT-GEE), make
vignette, and move the code relevant to permutation methods like
Score.test.stat.perm, Score.test.stat.zero.perm to CPP, complete the
roxygen skeleton.

## References:

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