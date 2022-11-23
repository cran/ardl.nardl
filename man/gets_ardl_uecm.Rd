\name{gets_ardl_uecm}
\alias{gets_ardl_uecm}
\title{
General-to-specific approach for the autoregressive distributed lag model
}
\description{
Adopt the general-to-specific approach to estimate the autoregressive distributed lag model
}
\usage{
gets_ardl_uecm(x, dep_var, expl_var, p_order = c(2), q_order = c(3),
              gets_pval = 0.1, case = 3, F_HC = FALSE, order_l = 5, 
              graph_save = FALSE)}
\arguments{
  \item{x}{
data frame
}
  \item{dep_var}{
A character vector. The dependent variable.
}
  \item{expl_var}{
Character vector. List of explanatory variable(s)
}
  \item{p_order}{
Integer. Maximum number of lags for 'dep_var' 
}
  \item{q_order}{
Integer. Maximum number of lags for 'expl_var'
}
  \item{gets_pval}{
Integer value between 0 and 1 needed for the general-to-specific approach. The default is 0.1 (10 percent significance level).  The chosen p-value is the criteria for determining non-significant repressors to be eliminated in a backward elimination path. The final parsimonious model is the best fit model based on the Schwarz information criteria
}
  \item{case}{
Positive integer 1 to 5. Default is 3
}
  \item{F_HC}{
Logical (default is FALSE). If TRUE, Heteroscedasticity-Consistent Covariance Matrix Estimation is applied to the model before when estimating F statistic
}
  \item{graph_save}{
Logical. If TRUE, display stability plot. Default is FALSE
}
  \item{order_l}{
Integer. order for the serial correlation, and heteroscedasticity test
}
}
\value{
\item{Parsimonious_ARDL_fit }{Return an estimated general-to-specific ARDL model}
\item{Parsimonious_ECM_fit }{Return an estimated general-to-specific error correction model}
\item{Summary_ecm_fit }{Return the summary of 'Parsimonious_ECM_fit'}
\item{Parsimonious_ECM_diagnostics_test }{Return the diagnostic test for 'Parsimonious_ECM_fit'.The diagnostic tests items are the Breusch-Godfrey test for higher-order serial correlation (BG_SC_lm_test). The Engle (1982) test for conditional heteroscedasticity (LM_ARCH_test). The test for non-normality is that of Jarque and Bera (1980). The RESET null hypothesis adopted implies - including the 2nd - degree terms improve the fit (over the model specified). Ljung and Box (1978) tests for autocorrelation in the residuals}
\item{cointegration }{Return the F statistic, the upper and lower critical values for PSS (2001) bounds test}
\item{Longrun_relation }{The estimated longrun relation from the error correction model}
}

\references{

Engle, R. F. (1982). Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflations. Econometrica (50) 987 - 1007

Ljung GM, Box GEP (1978). On a Measure of Lack of Fit in Time Series Models. Biometrika, 65(2), 297 - 303. https://doi.org/10.2307/2335207

Jarque C, Bera A (1980). Efficient Tests for Normality, Homoskedasticity, and Serial Independence. Economics Letters, 6(3), 255 - 259. https://doi.org/10.1016/0165-1765 (80) 90024-5

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships Journal of applied econometrics, 16(3), 289-326

}

\seealso{
\code{\link[ardl.nardl]{gets_nardl_uecm}}
\code{\link[ardl.nardl]{ardl_uecm}}
}
\examples{
  data(expectation)
  out <- gets_ardl_uecm(x = expectation,
                        dep_var = c('nq_inf_exp'), 
                        expl_var = c('food_inf','nethawkish'),
                        p_order = c(4), 
                        q_order = c(5,7),
                        gets_pval = 0.1, 
                        case = 4, 
                        graph_save = FALSE,
                        F_HC = FALSE,
                        order_l = 7)
  out
}
