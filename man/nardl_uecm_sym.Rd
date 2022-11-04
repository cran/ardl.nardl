\name{nardl_uecm_sym}
\alias{nardl_uecm_sym}
\title{
Short-run symmetry restrictions (SRSR) and the
long-run symmetry restrictions (LRSR).

}
\description{
~~
Performs one of Short-run symmetry restrictions or the
long-run symmetry restrictions.
}
\usage{
nardl_uecm_sym(x, assumption = c("SRSR"), decomp, control = NULL, c_q_order = NULL, 
p_order = c(3), q_order = c(4), dep_var, graph_save = FALSE, case = 3)
}
\arguments{
  \item{x}{
A dataframe object.  S3 dataframe object not applicable.
}
  \item{assumption}{
The Short-run symmetry restrictions (SRSR), the default, assume the transitory asymmetric effect is insignificant and testing the long-run asymmetric effect. Long-run symmetry restrictions (LRSR) imposes symmetric effect on the long-run coefficients of the positive and negative decomposed variable, and assuming the presence of transitory asymmetric effect.

}
  \item{decomp}{
A variable to be decomposed into positive (pos) and negative (neg)
}
  \item{control}{
A character vector with one element serving as the second variable single control variable. The default is NULL, but when control variable is specified, 'third' should be TRUE.
}
  \item{c_q_order}{
Lag order for the control. the default value is c(2) and only applicable when third = TRUE, and control is specified.
}
  \item{p_order}{
Number of lags for dependent variable (dep_var).
}
  \item{q_order}{
Short-run lags for each positive and negative decomposed variable. In case the LM test null hypothesis of no serial correlation is rejected, consider increasing the lags, in particular, the p_order or both p_order and q_order lags until the null hypothesis is accepted.
}
  \item{dep_var}{
Character vector which contain the dependent variable.
}
  \item{graph_save}{
Logical (default is FALSE). If TRUE, stability plot is displayed.}
  \item{case}{
A numerical vector of either 1, 2 3, 4 or 5. The default value is 3 - restricted constant (no intercept), no trend. See details for more information.
}
}
\details{
The Short-run symmetry restrictions (SRSR), the default may be useful in obtaining more accurate estimate when the sample is small (SYG, 2014). Besides, this may be necessary when the short-run asymmetric effect does not exist.
Long-run symmetry restrictions (LRSR) imposes symmetric effect on the positive and negative decomposed variable.
}
\value{
\item{NARDL_ECM_fit }{This display the NARDL ECM lm fit}
\item{UECM }{This display the summary results of the estimated unrestricted error correction models}
\item{cointegration }{This presents the Pesaran, Shin, and Smith (2001) cointegration test. The critical F and t statistics obtained if the sample is greater than 80 are asymptotic. If you adopt this segment, please consider citing Pesaran et al (2001) and Jordan and Philips (2020) dynamac package. If the sample is less than or equal to 80, for case 1, the F and t statistics adopted are large sample asymptotic critical value (CV) and no small-sample asymptotic critical value available. When case = 3 or 5, the critical value of F statistic is small sample asymptotic, and the adopted t statistics are large sample asymptotic CV. Consider including Narayan (2005) in your citation list. If the sample is greater than 80, for case 1, 3 and 5, large sample asymptotic CV for F and t statistics are adopted. In case 2 and 4, for sample <= 80, the F statistic are small sample asymptotic and no asymptotic CV for the t statistic. If sample > 80, for case 2 and 4, large sample asymptotic CV for F is adopted and no t-statistic CV.}
\item{Longrun_relation }{The estimated longrun relations}
\item{asymmetric test }{This contains the asymmetric test}
\item{diagnostics test}{This presents the Breusch-Godfrey test for higher-order serial correlation with null hypothesis of no serial correlation, the Autoregressive Conditional Heteroscedasticity (ARCH) LM-test with the null hypothesis of no ARCH effect and the Jarque-Bera normality test with null hypothesis of normality. This ARCH LM function is derived from Zaghdoudi (2018)}
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
\code{\link[ardl.nardl]{nardl_uecm}}
\code{\link[dynamac]{pssbounds}}
}
\examples{
data(syg_data)
out_srsr <- nardl_uecm_sym(x = syg_data,
             decomp  = 'ca_ip', 
             assumption = c('SRSR'),
             control =NULL,
             p_order =5,
             q_order =3,
             dep_var = 'ca_u',
             graph_save = FALSE, 
             case = 3)
out_srsr

out_lrsr <- nardl_uecm_sym(x = syg_data,
                           decomp  = 'ca_ip', 
                           assumption = c('LRSR'),
                           control =NULL,
                           p_order =5,
                           q_order =3,
                           dep_var = 'ca_u',
                           graph_save = FALSE, 
                           case = 3)
out_lrsr
}