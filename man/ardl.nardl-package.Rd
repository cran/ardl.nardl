\name{ardl.nardl-package}
\alias{ardl.nardl-package}
\alias{ardl.nardl}
\docType{package}
\title{
Linear and Nonlinear Autoregressive Distributed Lag Models: General-to-Specific Approach
}

\description{
Estimate the linear and nonlinear autoregressive distributed lag (ARDL & NARDL) models and the corresponding error correction models, and test for longrun and short-run asymmetric. The Pesaran, Shin & Smith (2001) Bounds test for level relationships is also provided with the aid of Jordan and Philips (2020) pssbounds function. In addition, the 'ardl.nardl' package also performs short-run and longrun symmetric restrictions available at Shin et al. (2014) and their corresponding tests.
}

\references{
Jordan S, Philips A (2020). _dynamac: Dynamic Simulation and Testing for Single-Equation
      ARDL Models_. R package version 0.1.11

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. Journal of applied econometrics, 16(3), 289-326.  https://doi.org/10.1002/jae.616
      
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling Asymmetric Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework. In: Sickles, R., Horrace, W. (eds) Festschrift in Honor of Peter Schmidt. Springer, New York, NY. https://doi.org/10.1007/978-1-4899-8008-3_9
}
\keyword{package}
\examples{
\dontrun{
data(fuel_price)
data(expectation)
 out1 <- gets_ardl_uecm(x = expectation,
                        dep_var = c('nq_inf_exp'), 
                        expl_var = c('food_inf','nethawkish'),
                        p_order = c(4), 
                        q_order = c(5,7),
                        gets_pval = 0.1, 
                        case = 4, 
                        graph_save = FALSE,
                        F_HC = FALSE,
                        order_l = 7)
 out1 

 out2 <- gets_nardl_uecm(x = expectation,
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
 out2

OUT3 <- auto_case_ardl(x = expectation,
               dep_var = 'n12m_inf_exp', 
               expl_var = c('food_inf',"hawkish","dovish"),
               p_order = 2,
               q_order = c(4,4,4),
               gets_pval = 0.05,
               graph_save = FALSE, 
               order_l = 7)
OUT3 
               
OUT4 <- nardl_auto_case(x = fuel_price,
                      decomp =  'wti',
                      dep_var = 'fpp',
                      control = 'bdc',
                      c_q_order = c(5), 
                      p_order = c(5),
                      q_order = c(6),
                      gets_pval = 0.1,
                      order_l = 4,
                      graph_save = FALSE)
OUT4

uecm_case3 <- ardl_uecm(x = fuel_price,
                            dep_var = c('fpp'),
                            expl_var = c('bdc', 'wti'),
                            p_order =c(6),
                            q_order =c(5,3),
                            graph_save = FALSE,
                            case = 3)
uecm_case3

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

data(syg_data)
out_srsr <- nardl_uecm_sym(x = syg_data,
             decomp  = 'ca_u', 
             assumption = c('SRSR'),
             control =NULL,
             p_order =5,
             q_order =3,
             dep_var = 'ca_ip',
             graph_save = FALSE, 
             case = 3)
out_srsr

out_lrsr <- nardl_uecm_sym(x = syg_data,
                           decomp  = 'ca_u', 
                           assumption = c('LRSR'),
                           control =NULL,
                           p_order =5,
                           q_order =3,
                           dep_var = 'ca_ip',
                           graph_save = FALSE, 
                           case = 3)
out_lrsr
}
}
