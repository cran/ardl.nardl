\name{ardl_uecm}
\alias{ardl_uecm}
\title{
Estimate the ARDL and ARDL Error Correction Model.
}
\description{
Estimate ARDL and ARDL error correction model, conduct the bounds test for cointegration 
}
\usage{
ardl_uecm(x, p_order = c(2), q_order, dep_var, graph_save = FALSE, expl_var, case = 3)
}
\arguments{
  \item{x}{
A dataframe object.  S3 dataframe object not applicable.
}
  \item{p_order}{
is the lag for the dependent variable (dep_var). Default value is c(2)
}
  \item{q_order}{
A numeric vector containing the lags for each explanatory variable. A vector of two lags (i.e c(3,4)) which corresponds to the expl_var (say c('x1','x2') where x1 has lag of 3 and x2 has lag of 4).
}
  \item{dep_var}{
Character vector containing one dependent variable
}
  \item{graph_save}{
Logical (default is FALSE). If TRUE, stability plot is displayed as well.
}
  \item{expl_var}{
Character vector containing the explanatory variable(s).
}
  \item{case}{
Vector of either 1, 2, 3, 4 or 5. The default value is 3 - restricted constant (no intercept), no trend
}
}
\value{
\item{ARDL_fit }{This display the ARDL lm fit}
\item{ARDL_ECM_fit }{This display the ARDL ECM lm fit}
\item{UECM }{This display the summary results of the estimated unrestricted error correction models}
\item{UECM }{This display the summary results of the estimated linear unrestricted error correction models}
\item{cointegration }{This presents the Pesaran, Shin, and Smith (2001) cointegration test. The critical F and t statistics obtained if the sample is greater than 80 are asymptotic. If you adopt this segment, please consider citing Pesaran et al (2001) and Jordan and Philips (2020) dynamac package. If the sample is less than or equal to 80, for case 1, the F and t statistics adopted are large sample asymptotic critical value (CV) and no small-sample asymptotic critical value available. When case = 3 or 5, the critical value of F statistic is small sample asymptotic, and the adopted t statistics are large sample asymptotic CV. Consider including Narayan (2005) in your citation list. If the sample is greater than 80, for case 1, 3 and 5, large sample asymptotic CV for F and t statistics are adopted. In case 2 and 4, for sample <= 80, the F statistic are small sample asymptotic and no asymptotic CV for the t statistic. If sample > 80, for case 2 and 4, large sample asymptotic CV for F is adopted and no t-statistic CV.}
\item{Longrun_relation }{The estimated longrun relations}
\item{diagnostics test}{This presents the Breusch-Godfrey test for higher-order serial correlation with null hypothesis of no serial correlation, the Autoregressive Conditional Heteroscedasticity (ARCH) LM-test with the null hypothesis of no ARCH effect and the Jarque-Bera normality test with null hypothesis of normality. This ARCH LM function is derived from Zaghdoudi (2018)}
}
\references{
Jordan S, Philips A (2020). _dynamac: Dynamic Simulation and Testing for Single-Equation /n
      ARDL Models_. R package version 0.1.11

Narayan, P. K. (2005). The saving and investment nexus for China: evidence from cointegration tests.
      Applied economics, 37(17), 1979-1990.

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level         relationships. Journal of applied econometrics, 16(3), 289-326.
      
Zaghdoudi, T. (2018). nardl: Nonlinear Cointegrating Autoregressive Distributed Lag Model_. R package version 0.1.5 
}
\seealso{
\code{\link[ardl.nardl]{nardl_uecm}}
\code{\link[ardl.nardl]{nardl_uecm_sym}}
}

\examples{
data(fuel_price)
uecm_case3 <- ardl_uecm(x = fuel_price,
                                   p_order =c(6),
                                   q_order =c(5,3),
                                   dep_var = c('fpp'),
                                   expl_var = c('bdc', 'wti'),
                                   graph_save = FALSE,
                                   case = 3)
uecm_case3
uecm_case3$cointegration
uecm_case3$Longrun_relation
uecm_case3$`diagnostics test`

uecm_case2 <- ardl_uecm(x = fuel_price,
                                  p_order =c(3),
                                  q_order =c(4),
                                  dep_var = c('fpp'),
                                  expl_var = c('wti'),
                                  graph_save = TRUE,
                                  case = 2)
uecm_case2
uecm_case2$cointegration
uecm_case2$Longrun_relation
uecm_case2$`diagnostics test`
}