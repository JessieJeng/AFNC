### AFNC
Adaptive False Negative Control (AFNC)

Z. John Daye and X. Jessie Jeng

AFNC: Performs the adaptive false negative control (AFNC), as descripted in Jeng et al. (2016).  The AFNC can provide informative analysis of very high-dimensional datasets with weak signals, such as in the analysis of large-scale sequencing studies at the single locus level, by including a large proportion of signals with high confidence.  The proportion of signals is estimated adaptively from the data.

This package requires a fortran 90 compiler.  Your build might fail if you don't have it.

#### Installation

Install AFNC from local source using

```r
install.packages("AFNC_1.0.tar.gz", repos=NULL, type="source")
```

Install AFNC from GitHubusing

```r
library(devtools)
install_github("zjdaye/AFNC")
```

#### Example use

The following example performs AFNC given a vector of p-values.

```r
set.seed(1)
cd = estimate.cd(d=length(p.value), alpha=0.05)
afnc = AFNC(p.value, alpha=0.05, beta=0.1, cd=cd)$afnc
selected = afnc$index # selected variables
selected.p.value = afnc$p.value # p-values of selected variables
```

The following example performs association test given a matrix of predictors and response.

```r
# Simulate response and predictors
set.seed(1); d = 10000; n = 2000
X = array(rnorm(n*d),c(n,d))
y = X[,1:50] %*% (1:50/10) + rnorm(n)

# Performs Wald's test
obj = association.test(X, y, method="Wald")
p.value = obj$p.value; test.stat = obj$test.stat
```

#### Licenses

The AFNC package as a whole is distributed under
[GPL-3 (GNU General Public License version 3)[http://www.gnu.org/licenses/gpl-3.0.en.html].
