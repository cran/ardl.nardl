\name{auto_case_ardl}
\alias{auto_case_ardl}
\title{
Obtain the best ARDL model specification and bounds test.
}
\description{
This function finds the best ARDL model specification and conduct bounds test by relying on the general to specific approach.
}
\usage{
auto_case_ardl(x, dep_var, expl_var, p_order, q_order, 
gets_pval = 0.05, order_l = 3, graph_save = FALSE)
}
\arguments{
  \item{x}{
Dataframe.
}
  \item{dep_var}{
A Character vector that contain the response variable.
}
  \item{expl_var}{
Character vector containing the list of explanatory variable(s).
}
  \item{p_order}{
An integer. Lag differenced adopted for the differenced response variable
}
  \item{q_order}{
  An integer. Lag differenced adopted for the differenced explanatory variable(s)
}
  \item{gets_pval}{
The p- value for adopted as criteria for eliminating non-significant variable in the course of obtaining the best model based on the Schwarz information criteria.
}
  \item{order_l}{
Integer. Needed for the autocorrelation and heteroscedasticity test 
}
  \item{graph_save}{
Logical. If TRUE, displays the stability plots
}
}
\details{
The procedure of the general-to-specific approach in obtaining the parsimonious model involves conducting the multi-path backwards elimination; tests both single and multiple hypothesis tests, diagnostics tests and goodness-of-fit measures. See page 5 of Sucarrat, (2021) for more details.

The value for gets_pval is influential the final model based on the multipath backward elimination. For more details on the general-to-specific approach, see the vignette of the 'gets' package.
}
\value{
\item{Parsimonious_ARDL_fit }{Return an estimated general-to-specific ARDL model}
\item{Parsimonious_ECM_fit }{Return an estimated general-to-specific error correction model}
\item{Summary_ecm_fit }{Return the summary of 'Parsimonious_ECM_fit'}
\item{Parsimonious_ECM_diagnostics_test }{Return the diagnostic test for 'Parsimonious_ECM_fit'.The diagnostic tests items are the Breusch-Godfrey test for higher-order serial correlation (BG_SC_lm_test). The Engle (1982) test for conditional heteroscedasticity (LM_ARCH_test). The test for non-normality is that of Jarque and Bera (1980). The RESET null hypothesis adopted implies - including the 2nd - degree terms improve the fit (over the model specified). Ljung and Box (1978) tests for autocorrelation in the residuals}
\item{cointegration }{Return the F statistic, the upper and lower critical values for PSS (2001) bounds test}
\item{Longrun_relation }{The estimated longrun relation from the error correction model}
}
\note{
Do not differenced the variables to be adopted in this function and all other functions for ARDL and NARDL estimation. The package inherently takes the difference and produced output with a prefix (D.) to the variable name and suffix the variable name with underscore (_) and the lag value.
}


\seealso{
\code{\link[gets]{gets}}
\code{\link[ardl.nardl]{gets_nardl_uecm}}
\code{\link[ardl.nardl]{ardl_uecm}}
\code{\link[ardl.nardl]{nardl_auto_case}}
}
\references{
Sucarrat, G. User-Specified General-to-Specific (GETS) and Indicator Saturation (ISAT) Methods. 28th September 2021. https://mirror.epn.edu.ec/CRAN/web/packages/gets/vignettes/user-defined-gets-and-isat.pdf
}
\examples{
data("expectation")
out_aut <- auto_case_ardl(x = expectation, 
                         dep_var = 'n12m_inf_exp', 
                         expl_var = c('food_inf',"hawkish","dovish"), 
                         p_order = 2, 
                         q_order = c(4,4,4), 
                         gets_pval = 0.05, 
                         graph_save = FALSE)

data("fuel_price")
out_aut <- auto_case_ardl(x = fuel_price, 
                         dep_var = 'fpp',
                         expl_var = c('bdc','wti'),
                         p_order = 2,
                         q_order = c(4,4),
                         gets_pval = 0.08,
                         graph_save = TRUE)
}
