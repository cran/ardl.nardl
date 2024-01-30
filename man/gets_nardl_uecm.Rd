\name{gets_nardl_uecm}
\alias{gets_nardl_uecm}
\title{
Parsimonious NARDL model
}
\description{
Adopt the general-to-specific approach to estimate the autoregressive distributed lag model
}
\usage{
gets_nardl_uecm(x, decomp, dep_var, control = NULL, d = Inf, c_q_order = c(2), 
                p_order = c(3), q_order = c(4), gets_pval = 0.1, order_l = 4, 
                F_HC = FALSE, graph_save = FALSE, case = 3)}
\arguments{
  \item{x}{
data frame
}
  \item{decomp}{
A character vector. The variable decomposed into positive (pos) and negative (neg) change.
}
  \item{dep_var}{
  A character vector. The dependent variable 
}
  \item{control}{
A character vector. Default is NULL. The second dependent variable.
}
  \item{d}{
An integer or character vector. The threshold value (d) is adopted when computing the partial sum. The default value Inf is the implicit threshold for the partial sum whenever the base 'cumsum' function is adopted. The value of d can be 'mean' or 0 and any other integer specified as the threshold value for the partial sum.
}
  \item{c_q_order}{
Integer. Maximum number of lags for 'control'
}
  \item{p_order}{
Integer. Maximum number of lags for 'dep_var' 
}
  \item{q_order}{
Integer. Maximum number of lags for level and differenced 'decomp'
}
  \item{gets_pval}{
Integer value between 0 and 1 needed for the general-to-specific approach. The default is 0.1 (10 percent significance level).  The chosen p-value is the criteria for determining non-significant repressors to be eliminated in a backward elimination path. The final parsimonious model is the best fit model based on the Schwarz information criteria. 
}
  \item{order_l}{
Integer. order for the serial correlation, and heteroscedasticity test
}
  \item{F_HC}{
Logical (default is FALSE). If TRUE, heteroscedasticity-Consistent Covariance Matrix Estimation is applied to the model before when estimating F statistic
}
  \item{graph_save}{
Logical. If TRUE, display stability plot. Default is FALSE.
}
  \item{case}{
Positive integer 1 to 5. Default is 3
}
}
\value{
\item{Parsimonious_NARDL_fit }{Return an estimated general-to-specific NARDL model.}
\item{Parsimonious_ECM_fit }{Return an estimated general-to-specific error correction model.}
\item{Summary_uecm_fit }{Return the summary of 'Parsimonious_ECM_fit'}
\item{ecm_diagnostics_test }{Return the diagnostic test for the 'Parsimonious_ECM_fit'. The diagnostic tests indicate the Breusch-Godfrey test for higher-order serial correlation (BG_SC_lm_test). The Engle (1982) test for conditional heteroscedasticity (LM_ARCH_test). The test for non-normality is that of Jarque and Bera (1980). The RESET null hypothesis adopted implies - including the 2nd - degree terms improve the fit (over the model specified).}
\item{longrun_asym }{Return the estimated longrun asymmetric test}
\item{Shortrun_asym }{Return the estimated short-run asymmetric test.If one of the decomposed variable does not appear among the shortrun differenced variables of the parsimonious model, The value returned is a wald test of whether the sum of the coefficients of the remaining decomposed variable included does not have any significant effect on the best model }
\item{cointegration }{Return the F statistic, the upper and lower critical values for PSS (2001) bounds test. Please, disregard the tstat on the cointegration test.}
\item{Longrun_relation }{The longrun relation}
}
\references{
Engle, R. F. (1982). Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflations. Econometrica 50: 987 - 1007.

Jarque C, Bera A (1980). Efficient Tests for Normality, Homoskedasticity, and Serial Independence. Economics Letters, 6(3), 255 - 259. https://doi.org/10.1016/0165-1765(80) 90024-5.

Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling Asymmetric Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework. In: Sickles, R., Horrace, W. (eds) Festschrift in Honor of Peter Schmidt. Springer, New York, NY. https://doi.org/10.1007/978-1-4899-8008-3_9

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. Journal of applied econometrics, 16(3), 289-326.
}
\seealso{
\code{\link[ardl.nardl]{gets_ardl_uecm}}
\code{\link[ardl.nardl]{ardl_uecm}}
}
\note{
The decomposed variable should display both positive and negative change, preferably on a balanced scale. However, when a variable display only positive change and no negative change, vice versa, such variable should not be adopted (i.e decomposed).
}
\examples{
\dontrun{
  data(expectation)
  out <- gets_nardl_uecm(x = expectation,
                         decomp = 'food_inf',
                         dep_var = 'nq_inf_exp',
                         control = 'nethawkish',
                         c_q_order = c(3),
                         p_order = c(3),
                         q_order = c(3),
                         gets_pval = 0.1,
                         graph_save = FALSE,
                         case = 5, 
                         F_HC = FALSE)
  out
}
}
