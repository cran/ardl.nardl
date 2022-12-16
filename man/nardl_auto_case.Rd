\name{nardl_auto_case}
\alias{nardl_auto_case}
\title{
Obtain the best NARDL model specification and bounds test.
}
\description{
This function finds the best NARDL model specification and conduct bounds test by relying on the general to specific approach.
}
\usage{
nardl_auto_case(x, decomp, dep_var, control = NULL, c_q_order = c(2), 
                p_order = c(3), q_order, gets_pval = 0.1, 
                order_l = order_l, graph_save = FALSE)
}
\arguments{
  \item{x}{
Dataframe
}
  \item{decomp}{
A character vector. The variable to be decomposed to positive (pos) and negative (neg) variable.
}
  \item{dep_var}{
  A character vector. The dependent variable 
}
  \item{control}{
A character vector. Default is NULL. The second dependent variable.
}
  \item{c_q_order}{
Integer. Maximum number of lags for 'control'
}
  \item{p_order}{
An integer. Lag differenced adopted for the differenced response variable
}
  \item{q_order}{
  An integer. Lag differenced adopted for the differenced explanatory variable(s)
}
  \item{gets_pval}{
The p- value which served as criteria for eliminating non-significant variable in the course of obtaining the best model based on the Schwarz information criteria.
}
  \item{order_l}{
Integer. Needed for the autocorrelation and heteroscedasticity test 
}
  \item{graph_save}{
Logical. If TRUE, displays the stability plots
}
}
\details{
The procedure of
the general-to-specific approach in obtaining the parsimonious model involves conducting the multi-path backwards elimination; tests both single and multiple hypothesis tests, diagnostics tests and goodness-of-fit measures. See page 5 - 6 of Sucarrat (2021) for more details.

The value for gets_pval is influential the final model based on the multipath backward elimination. For more details on the general-to-specific approach, see the vignette of the 'gets' package.
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
\note{
Do not differenced the variables to be adopted in this function and all other functions for NARDL and ARDL estimation. The package inherently takes the difference and produced output with a prefix (D.) to the variable name and suffix the variable name with underscore (_) and the lag value. 
}

\seealso{
\code{\link[gets]{gets}}
\code{\link[ardl.nardl]{gets_nardl_uecm}}
\code{\link[ardl.nardl]{ardl_uecm}}
\code{\link[ardl.nardl]{auto_case_ardl}}
}
\references{
Sucarrat, G. User-Specified General-to-Specific (GETS) and Indicator Saturation (ISAT) Methods. 28th September 2021. https://mirror.epn.edu.ec/CRAN/web/packages/gets/vignettes/user-defined-gets-and-isat.pdf
}
\examples{
\dontrun{
 data("fuel_price")
 out1 <- nardl_auto_case(x = fuel_price,
                decomp =  'wti',
                dep_var = 'fpp',
                control = 'bdc',
                c_q_order = c(5), 
                p_order = c(5),
                q_order = c(6),
                gets_pval = 0.1,
                order_l = 4,
                graph_save = FALSE)
 out1

 out2 <- nardl_auto_case(x = fuel_price,
                decomp =  'wti',
                dep_var = 'fpp',
                control = NULL,
                c_q_order = c(4), 
                p_order = c(5),
                q_order = c(6),
                gets_pval = 0.02,
                order_l = 4,
                graph_save = FALSE)
 out2
}
}
