\name{nardl_uecm}
\alias{nardl_uecm}
\title{
Estimate the nonlinear ARDL (NARDL) Error Correction Model
}
\description{
Function to estimate Shin, Yu, and Greenwood-Nimmo (2014) nonlinear ARDL (NARDL), the NARDL ECM, and conduct the bounds test for cointegration.
}
\usage{
nardl_uecm(x, decomp, d = Inf, control = NULL, c_q_order = c(2), 
p_order = c(3), q_order = c(4), dep_var, order_l = 4, graph_save = FALSE, case = 3)
}
\arguments{
  \item{x}{
A dataframe object.  S3 dataframe object not applicable.
}
  \item{decomp}{
A character vector. The variable to be decomposed into positive (pos) and negative (neg) change.
}
  \item{d}{
An integer or character vector. The default value Inf for the threshold value (d) is adopted when computing the partial sum. Inf is the implicit threshold for the partial sum whenever the base 'cumsum' function is adopted. The value of d can be 'mean' or 0 and any other integer specified as the threshold value for the partial sum.
}
\item{control}{
A character vector with one element serving as the second variable in the right hand side of the equation. The default is NULL.
}
\item{c_q_order}{
Lag order for 'control'. the default value is c(2) and only applicable when control is specified.
}
\item{p_order}{
Positive interger. Number of lags for dependent variable (dep_var).
}
\item{q_order}{
Positive interger. Short-run lags for the decomposed variables. In case the LM test null hypothesis of no serial correlation is rejected, consider increasing the lags, in particular, the p_order or both p_order and q_order lags until the null hypothesis is accepted.
}
\item{dep_var}{
Character vector. The dependent variable.
}
\item{order_l}{
Positive integer. lag order for the serial correlation, and heteroscedasticity test
}
\item{graph_save}{
Logical (default is FALSE). If TRUE, stability plot is displayed.}
\item{case}{
A numerical vector that can take a value between 1 to 5. The default value is case 3. See details for more information.
}
}
\value{
\item{NARDL_fit }{This display the NARDL lm fit}
\item{NARDL_ECM_fit }{This display the NARDL ECM lm fit}
\item{UECM }{This display the summary results of the estimated unrestricted error correction models}
\item{cointegration }{This presents the Pesaran, Shin, and Smith (2001) cointegration test. The critical F and t statistics obtained if the sample is greater than 80 are asymptotic. If you adopt this segment, please consider citing Pesaran et al (2001) and Jordan and Philips (2020) dynamac package. If the sample is less than or equal to 80, for case 1, the F and t statistics adopted are large sample asymptotic critical value (CV) and no small-sample asymptotic critical value available. When case = 3 or 5, the critical value of F statistic is small sample asymptotic, and the adopted t statistics are large sample asymptotic CV. Consider including Narayan (2005) in your citation list. If the sample is greater than 80, for case 1, 3 and 5, large sample asymptotic CV for F and t statistics are adopted. In case 2 and 4, for sample <= 80, the F statistic are small sample asymptotic and no asymptotic CV for the t statistic. If sample > 80, for case 2 and 4, large sample asymptotic CV for F is adopted and no t-statistic CV.}
\item{Longrun_relation }{The estimated longrun relations}
\item{Longrun_asymmetric_test }{This contains the asymmetric test, with the null hypothesis of longrun coefficients of the positive decomposed variable equals to the coefficient of negative decomposed variable.}
\item{Shortrun_asymmetric_test }{This contains the asymmetric test, with the null hypothesis the sum of the coefficients of all the short run positive decomposed variable equals to the sum of the coefficients of the short-run negative decomposed variable.}
\item{diagnostics test}{
The diagnostic tests indicate the Breusch-Godfrey test for higher-order serial correlation (BG_SC_lm_test). The Engle (1982) test for conditional heteroscedasticity (LM_ARCH_test). The test for non-normality is that of Jarque and Bera (1980). The RESET null hypothesis adopted implies - including the 2nd - degree terms improve the fit (over the model specified). This ARCH LM function is derived from Zaghdoudi (2018)}
}
\references{
Jordan S, Philips A (2020). _dynamac: Dynamic Simulation and Testing for Single-Equation
      ARDL Models_. R package version 0.1.11

Narayan, P. K. (2005). The saving and investment nexus for China: evidence from cointegration tests.
      Applied economics, 37(17), 1979-1990.

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level         relationships. Journal of applied econometrics, 16(3), 289-326.
      
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling Asymmetric Cointegration and Dynamic              Multipliers in a Nonlinear ARDL Framework. In: Sickles, R., Horrace, W. (eds) Festschrift in          Honor of Peter Schmidt. Springer, New York, NY. https://doi.org/10.1007/978-1-4899-8008-3_9

Zaghdoudi, T. (2018). nardl: Nonlinear Cointegrating Autoregressive Distributed Lag Model_. R package version 0.1.5 

}
\seealso{
\code{\link[ardl.nardl]{ardl_uecm}}
\code{\link[ardl.nardl]{nardl_uecm_sym}}
}
\note{
The decomposed variable should display both positive and negative trend, preferably on a balanced scale. However, when a variable display only positive change and no negative change, vice versa, such variable should not be adopted (i.e decomposed).
}
\examples{
data(fuel_price)
output_n1_case5 <- nardl_uecm(x = fuel_price,
                                 decomp  = c('bdc'), 
                                 control =c('wti'),
                                 c_q_order = c(2), 
                                 p_order = c(3), 
                                 q_order = c(5), 
                                 dep_var = c('fpp'), 
                                 graph_save = FALSE, 
                                 case = 5)
output_n1_case5
output_n1_case5$UECM

output_n2_case3 <- nardl_uecm(x = fuel_price, 
                                 decomp  = c('bdc'), 
                                 control =NULL,
                                 c_q_order = NULL, 
                                 p_order = c(3), 
                                 q_order = c(5), 
                                 dep_var = c('fpp'), 
                                 graph_save = FALSE, 
                                 case = 3)
output_n2_case3
}