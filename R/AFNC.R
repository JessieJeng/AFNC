######################################################
# AFNC
#' Performs AFNC
#'
#' AFNC takes a vector of p-values and determines the variables (SNPs) selected by AFNC.  The AFNC can retain a high proportion of signals by discarding variants from the Noise region using adaptive false negative control (Jeng et al. 2016). Proportion of signals is estimated from the data.
#'
#' @param p.value d-vector of p-values.
#' @param alpha level of family-wise error rate for false positive control using the Bonferroni.
#' @param beta level of signal missing rate for false negative control using the AFNC.
#' @param cd \eqn{c_d} for controlling Type I error rate under the global null, must be pre-computed using {estimate.cd}.
#' @param c0 maximum number of the most significant p-values considered.  If {NULL}, {c0} will be set to \eqn{\min( \max(5000,0.1*d), d )}, such that at least 10\% of the p-values are considered.
#'
#' @return
#'  \item{bonferroni}{Two vectors ({index}, {p.value}).  {index} is a vector of indices of variables (SNPs) selected by the Bonferroni, and {p.value} is the corresponding p-values.}
#'  \item{afnc}{Two vectors ({index}, {p.value}).  {index} is a vector of indices of variables (SNPs) selected by the AFNC, and {p.value} is the corresponding p-values.}
#'  \item{t.alpha}{Two values ({rank}, {p.value}).  {rank} specifies the rank by which all variants ranked at or before \eqn{t_\alpha} are retained for controlling false positives by Bonferroni.  {p.value} is the corresponding p-value, such that all variants less than or equal to {p.value} are retained.}
#'  \item{T.fn}{Two values ({rank}, {p.value}). {rank} specifies the rank by which all variants ranked at or before \eqn{t_\alpha} are retained for adaptive false negative control (AFNC).  {p.value} is the corresponding p-value, such that all variants less than or equal to {p.value} are retained.}
#'  \item{signal.proportion}{The estimated signal proportion \eqn{\hat{\pi}}.}
#'  \item{number.signals}{The estimated number of signals \eqn{\hat{s} = \hat{\pi} * d}.}
#'
#' @details
#' The algorithm implemented in this function is as follows.  (See Jeng et al. (2016) for further details.)
#' \enumerate{
#'   \item The p-values are ordered at decreasing significance.
#'   \item The signal proportion estimator \eqn{\hat{\pi}} and estimated number of signals \eqn{\hat{s} = \hat{\pi} * d} are obtained using {estimate.signal.proportion}.
#'   \item Two cutoff positions, \eqn{t_\alpha} and \eqn{T_fn}, are determined to separate the Signal, Indistinguishable, and Noise regions. (See Figure 1 of Jeng et al. (2016) for illustration of the Signal, Indistinguishable, and Noise regions of inference.)
#'   \item Finally, variables (SNPs) with ordered p-values ranked at or before \eqn{t_\alpha} are selected by Bonferroni for family-wise false positive control.  Variables (SNPs) with ordered p-values ranked at or before \eqn{T_fn} are selected by the AFNC procedure for adaptive false negative control.
#' }
#'
#' @references
#'  Jeng, X.J., Daye, Z.J., Lu, W., and Tzeng, J.Y. (2016) Rare Variants Association Analysis in Large-Scale Sequencing Studies at the Single Locus Level.
#'
#' @examples
#' # Estimate signal proportions
#' set.seed(1)
#' cd = estimate.cd(d=length(p.value), alpha=0.05)
#' afnc = AFNC(p.value, alpha=0.05, beta=0.1, cd=cd)$afnc
#' selected = afnc$index # selected variables
#' selected.p.value = afnc$p.value # p-values of selected variables
#'
#' @export
######################################################
AFNC = function(p.value, alpha=0.05, beta=0.1, cd, c0=NULL) {
  d = length(p.value); pnames = names(p.value)

  # Maximum number of signals to consider
  if (is.null(c0)) { c0=max(5000,.1*d) }
  c0 = min(c0,d);

  # Estimates cd
  if (is.null(cd)){ cd  = estimate.cd(d, M=M, alpha) }

  # Use original ranking as names to p.values
  if (is.null(pnames)) { pnames = paste("V", 1:d, sep=""); names(p.value) = pnames }

  #-----------------------------
  # 1) Sort p-values
  sort.obj = sort(p.value, method="quick", index.return = TRUE)
  sorted.pvalue = sort.obj$x; sorted.index = sort.obj$ix; names(sorted.index) = pnames[sorted.index]

  #-----------------------------
  # 2) Estimate signal proportion
  signal.obj = estimate.signal.proportion(sorted.pvalue, cd, c0, is.sorted=TRUE)
  signal.prop = signal.obj$signal.proportion; number.signals = signal.obj$number.signals

  #-----------------------------
  # 3) Determine cutoff positions
  # Determine 1st cut for family-wise Type I error rate control
  cut1 = max( which(sorted.pvalue < alpha/((1-signal.prop)*d)), 0 )

  # Determine 2nd cut for false negative control
  p_j = sorted.pvalue[floor(signal.prop*d+1) : d]
  ind=1:d; numNul=length(p_j)
  quant = qbeta(beta, ind[1: numNul], numNul-(ind[1: numNul])+1)
  if (round(signal.prop*d) == 0){  # If no signals
    cut2 = 0; cut1 = 0
  } else{
    if (round(signal.prop*d) <= cut1) {  # cutoff positions coincides
      cut2 = cut1
    } else{
      jhat =min(which(p_j <= quant), numNul);
      cut2 = floor(signal.prop*d+jhat)
    }
  }

  #-----------------------------
  # 4) Select variables
  # Return values
  cutoff_p.value = c(0,0)
  bonferroni = NULL; afnc = NULL; fdr = NULL
  bonferroni.pval = NULL; afnc.pval = NULL; fdr.pval = NULL

  # Bonferroni
  if(cut1>0) {
    bonferroni = sorted.index[1:cut1]
    bonferroni.pval = sorted.pvalue[1:cut1]
    cutoff_p.value[1] = sorted.pvalue[cut1]
  }

  # afnc
  if(cut2>0) {
    afnc = sorted.index[1:cut2]
    afnc.pval = sorted.pvalue[1:cut2]
    cutoff_p.value[2] = sorted.pvalue[cut2]
  }

  return(list(bonferroni=list(index=bonferroni, p.value=bonferroni.pval),
              afnc=list(index=afnc, p.value=afnc.pval),
              t.alpha=list(rank=cut1, p.value=cutoff_p.value[1]),
              T.fn=list(rank=cut2, p.value=cutoff_p.value[2]),
              signal.proportion=signal.prop, number.signals=number.signals))
}

#====================================================================
# Association test
#' Performs association test
#'
#' {association.test} performs association tests at each predictor and returns both the test statistics and p-values.  The Wald test can be employed for quantitative traits.  The score or Lagrange multiplier test can be employed for both quantitative and qualitative traits. (See `Details'.)
#'
#' @param X n-by-d matrix of predictors (genotypes) for n samples and d variables (SNPs).
#' @param y Response (phenotype / (disease) trait).
#' @param method Association test to be used. See `Details'.
#' @param trait Quantitative or qualitative trait.  If it is {NULL}, qualitative is automatically used if there are 2 unique response values, and quantitative otherwise.
#' @return
#'    \item{test.stat}{d-vector of test statistics at each SNP.}
#'    \item{p.value}{d-vector of p-values at each SNP.}
#' @details
#'  Method "{score}" uses a Lagrange multiplier or score test, computed under the null model.
#'  Method "{Wald}" uses a Wald's test or t-test for quantitative traits, computed under the unrestricted model.  Wald's test is not recommended for qualitative traits, due to potential inefficiency of the Wald's test under logistic regression (Hauck and Donner 1977).
#'
#'  Trait "{qualitative}" is based on logistic regression.
#'  Trait "{quantitative}" is based on linear regression.  See Kim et al. (2014).
#'
#' @references
#'  Kim, S., Pan W., and Shen, X. (2014) Penalized regression approaches to testing for quantitative trait-rare variant association.  {Frontiers in Genetics}, 5:121.
#'
#'  Hauck, W.W. and Donner, A. (1977) Wald's test as applied to hypotheses in logit analysis. {J Am Stat Assoc}, 72, 851-853.
#' @examples
#' # Load "AFNC" library and example data.
#' library("AFNC")
#' data(example_data)
#'
#' # Simulate response and predictors
#' set.seed(1); d = 10000; n = 2000
#' X = array(rnorm(n*d),c(n,d))
#' y = X[,1:50] %*% (1:50/10) + rnorm(n)
#'
#' # Performs Wald's test
#' obj = association.test(X, y, method="Wald")
#' p.value = obj$p.value; test.stat = obj$test.stat
#'
#' @useDynLib AFNC
#' @export
#====================================================================
association.test <- function(X, y, method="score", trait=NULL) {
  n = dim(X)[1]; d = dim(X)[2]
  test.stat = rep(NA, d); p.value = rep(NA, d)
  snpNames = colnames(X)

  # Determine if qualitative
  if (is.null(trait)) { trait = ifelse( length(unique(y))==2, "qualitative", "quantitative" ) }

  # Recode into 0, 1 with the larger response value as 1.
  if (trait=="qualitative") { if (any(c(0,1)!=sort(unique(y)))) {ytmp=rep(0,length(y)); ytmp[y==max(y)]=1; y=ytmp } }

  # Wald's test
  if (method=="Wald") {
    # Quantitative
    if (trait=="quantitative") {
      #for (j in 1:d) { r = cor(X[,j],y); test.stat[j] = r*sqrt(n-2)/sqrt(1-r^2) }
      #test.stat.old = test.stat
      test.stat = .Fortran("FWald_quantitative", X=as.double(X), y=as.double(y),
               n=as.integer(n), d=as.integer(d), test_stat=numeric(d), package="AFNC")$test_stat

      #if ( any(test.stat!=test.stat.old) ) { save.image(file="tmp.RData"); stop("test.stat different") }
    } else { # Case-control
      stop("Wald's test not recommended for qualitative trait.")
    }
    p.value = 2*pt(abs(test.stat), df=n-2, lower.tail=F)
  }

  # Score test
  if (method=="score") {
    # Quantitative
    if (trait=="quantitative") {
      test.stat = .Fortran("Fscore_quantitative", X=as.double(X), y=as.double(y),
                           n=as.integer(n), d=as.integer(d), test_stat=numeric(d), package="AFNC")$test_stat
    } else { # Case-control
      test.stat = .Fortran("Fscore_qualitative", X=as.double(X), y=as.double(y),
                           n=as.integer(n), d=as.integer(d), test_stat=numeric(d), package="AFNC")$test_stat
    }
    test.stat[ is.na(test.stat) | is.infinite(test.stat) | is.nan(test.stat) ] = 0
    p.value = pchisq(test.stat, df=1, lower.tail=F)
  }

  names(test.stat) = snpNames; names(p.value) = snpNames
  return(list(p.value=p.value, test.stat=test.stat))
}

#====================================================================
# estimate.cd
#' Estimates \eqn{c_d} for controlling Type I error rate under the global null.
#'
#' {estimate.cd} empirically estimates \eqn{c_d} to control the Type I error rate under the global null hypothesis that no signals exist.  In this algorithm, \eqn{M} number of Monte Carlo samples, drawn under the global null hypothesis of no signals, are used to estimate \eqn{c_d}.  The estimate of \eqn{c_d} is pre-computed in {estimate.cd} and passed to {AFNC}.
#'
#' @param d number of variables (SNPs).
#' @param M number of Monte Carlo samples. The larger the \eqn{M} the more accurate the estimate.
#' @param alpha significance level for false positive control.
#'
#' @return Estimate of \eqn{c_d} for controlling Type I error rate under the global null, to be used as an input in {AFNC}.
#'
#' @examples
#' # Load "AFNC" library and example data.
#' library("AFNC")
#' data(example_data)
#'
#' set.seed(1)  # Set seed
#' cd = estimate.cd(d=length(p.value), M=10000, alpha=0.05)  # estimate c_d
#'
#' @export
#====================================================================
estimate.cd <- function(d, M=10000, alpha=0.05) {
  Vn = rep(0, M)
  idxProp = (1:d) / d
  for(i in 1 : M){
    noise_p <- pnorm(rnorm(d), lower.tail=F)
    sorted.pvalue = sort(noise_p, method="quick", index.return = FALSE)
    Vn[i] = max( ( idxProp - sorted.pvalue ) / sqrt( sorted.pvalue * (1 - sorted.pvalue) ) )
  }
  cd <- quantile(Vn, 1 - alpha)
  return(cd);
}

######################################################
# estimate.signal.proportion
#' Estimates signal proportions
#'
#' {estimate.signal.proportion} estimates the signal proportion \eqn{\pi} using a modified estimator based on Meinshausen and Rice (2006). The estimator does not depend on statistical normality assumptions and can adapt to \eqn{\pi} small.  (See Jeng et al. (2016) for details.)
#'
#' @param p.value d-vector of p-values
#' @param cd \eqn{c_d} for global Type I error control
#' @param c0 maximum number of the most significant p-values considered.  If {NULL}, {c0} will be set to \eqn{\min( \max(5000,0.1*d), d )}, such that at least 10\% of the p-values are considered.
#' @param is.sorted if p.value is already sorted.  If {TRUE}, p.value will not be re-sorted.
#'
#' @return
#'  \item{signal.proportion}{The estimated signal proportion \eqn{\hat{\pi}}.}
#'  \item{number.signals}{The estimated number of signals \eqn{\hat{s} = \hat{\pi} * d}.}
#'
#' @references
#'  Jeng, X.J., Daye, Z.J., Lu, W., and Tzeng, J.Y. (2016) Rare Variants Association Analysis in Large-Scale Sequencing Studies at the Single Locus Level.
#'
#'  Meinshausen M. and Rice, J. (2006) Estimating the proportion of false null hypotheses among a large number of independent tested hypotheses.  {Ann. Statist.}, 34:373-393.
#'
#' @examples
#' # Load "AFNC" library and example data.
#' library("AFNC")
#' data(example_data)
#'
#' # Estimate signal proportions
#' set.seed(1)
#' cd = estimate.cd(d=length(p.value), M=10000, alpha=0.05) # Estimate cd
#' obj = estimate.signal.proportion(p.value, cd=cd)
#' signal.proportion = obj$signal.proportion  # Estimated signal proportion
#' number.signals = obj$number.signals  # Estimated number of signals
#'
#' @export
######################################################
estimate.signal.proportion = function(p.value, cd, c0=NULL, is.sorted=FALSE) {
  d = length(p.value)

  # Maximum number of signals to consider
  if (is.null(c0)) { c0=max(5000,.1*d) }
  c0 = min(c0,d);

  # Sort p-values by decreasing significance
  if (is.sorted) { sorted.pvalue=p.value } else { sorted.pvalue=sort(p.value, method="quick", index.return = FALSE) }

  # Estimates the signal proportion \hat{\pi}
  W = ( (1:d)/d - sorted.pvalue - cd*sqrt(sorted.pvalue*(1-sorted.pvalue)))/(1-sorted.pvalue)
  signal.prop = max( W[2: c0], 0 )

  return( list(signal.proportion=signal.prop, number.signals=signal.prop*d) )
}

######################################################
#' Example data
#'
#' This dataset contains 10,000 p-values.  These data are used to illustrate the functions {AFNC} and {estimate.signal.proportion}.
#'
#' @docType data
#' @keywords datasets
#' @name example_data
#' @usage data(example_data)
#' @format A data set with a vector p-values.
#'
#' @examples
#' # Load "AFNC" library and example data.
#' library("AFNC")
#' data(example_data)
#'
#'
######################################################
NULL

