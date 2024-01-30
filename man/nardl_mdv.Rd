\name{nardl_mdv}
\alias{nardl_mdv}
\title{
A NARDL model with two decomposed variables
}
\description{
Estimate a NARDL model, having mulitple (two) decomposed variables
}
\usage{
nardl_mdv(x, dep_var, decomp1, decomp2, thresh1 = Inf, thresh2 = Inf, 
gets = TRUE, gets_pval = 0.1, case = NULL, conservative = FALSE, p_order = c(3), 
q_order1 = c(5), q_order2 = c(5), order_l = 4, graph_save = FALSE)
}
\arguments{
  \item{x}{
A dataframe
}
  \item{dep_var}{
The dependent variable
}
  \item{decomp1}{
Initial variable to be decomposed into postive and negative.  
}
  \item{decomp2}{
The second variable decomposed into postive and negative.  
}
  \item{thresh1}{
An integer or character vector. The threshold value (thresh1) for the first variable is adopted when computing the partial sum. Inf is the implicit threshold value for the partial sum whenever the base 'cumsum' function is adopted. The value of the threshold can be 'mean' or 0 or any other integer stated for the partial sum.
}
  \item{thresh2}{
An integer or character vector. The threshold value (thresh2) for the second decomposed variable.
}
  \item{gets}{
Logical. General-to-specific (GETS) approach. Default is TRUE which indicate adopting the GETS
}
  \item{gets_pval}{
The p-value adopted when gets is set as TRUE
}
  \item{case}{
An integer which can take either of 1, 2, 3, 4 or 5. The assumption on Bounds test procedure. Default is NULL. When it is set as an integer, gets should be set as FALSE.
}
  \item{conservative}{
Logical. Default is FALSE. When TRUE, the decomposed variables are assumed to be k = 2. When FALSE, k = 4.
}
  \item{p_order}{
An integer. Take the number of lags applicable to the dependent variable
}
  \item{q_order1}{
An integer. Take maximum number of lags applicable to the first variable
}
  \item{q_order2}{
An integer. The maximum number of lags applicable to the second variable}
  \item{order_l}{
An integer. Used in the diagnostics test
}
  \item{graph_save}{
Logical. Default is FALSE. When TRUE, return the stability plots of the model
}
}
\details{
Return a list containing 
}
\value{
\item{NARDL_fit }{NARDL model}
\item{ECM_fit }{ NARDL-ECM}
\item{Summary_uecm_fit }{Summary of ECM_fit}
\item{ecm_diagnostics_test }{Diagnostic tests}
\item{longrun_asym }{longrun asymmetric test}
\item{Shortrun_asym }{Shortrun asymmetric test}
\item{cointegration }{PSS bounds test}
\item{Longrun_relation }{Longrun relationship}
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
\code{\link[ardl.nardl]{nardl_uecm}}
}
\note{
The decomposed variable should display both positive and negative change, preferably on a balanced scale. However, when a variable display only positive change and no negative change, vice versa, such variable should not be adopted (i.e decomposed).
}
\examples{
\dontrun{
data(expectation)
nardl_mdv(x = expectation,
          dep_var = 'nq_inf_exp',
          decomp1 = 'food_inf',
          decomp2 = 'nethawkish',
          p_order = c(7),
          q_order1 = c(4),
          q_order2 = c(6),
          gets_pval = 0.1,
          conservative = FALSE,
          gets = FALSE,
          case = 5,
          order_l = 3,
          graph_save = FALSE)
}}
