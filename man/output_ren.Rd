\name{output_ren}
\alias{output_ren}
\title{
Rename the coefficient of the nardl_uecm summary output
}
\description{
As the case might be, their could be need to rename the rownames in nardl_uecm summary output
}
\usage{
output_ren(x, listn = listn, D.patern, D.repl, l.patern, l.repl)
}
\arguments{
  \item{x}{
A summary model output, preferably, nardl_uecm object
}
  \item{listn}{
A character vector

}
  \item{D.patern}{
The difference and lag difference variables contained the model
}
  \item{D.repl}{
Character vector containing the new names to replace the D.patern
}
  \item{l.patern}{
Character vector of variable appearing at level
}
  \item{l.repl}{
Character vector containing the names to replace l.patern
}
}
\value{
\item{x }{nardl_uecm or nardl_uecm_sym object}
}

\examples{
data(syg_data)
un <- c('ca_u','jp_u','us_u')
ip <- c('ca_ip','jp_ip','us_ip')
cnames <- c('Canada','Japan','United States')
case <- c(5,5,5)
p_lags <- c(3,4,5)
q_lags <- c(3,4,5)
########################################
out_nardl <- lapply(1:length(un), function(i)
  nardl_uecm(x = syg_data,
             decomp  = un[[i]], 
             control =NULL,
             c_q_order = c(2), 
             p_order =p_lags[[i]],
             q_order =q_lags[[i]],
             dep_var = ip[[i]],
             graph_save = FALSE, 
             case = 3))
names(out_nardl) <- cnames
out_nardl$Canada$UECM

output <- output_ren(x=out_nardl,listn = un,
           D.patern="D.ca_ip|D.jp_ip|D.us_ip",D.repl="D.output",
           l.patern ="ca_ip|jp_ip|us_ip",l.repl="output")

output <- output_ren(x=output,listn = un,
             D.patern="D.ca_u|D.jp_u|D.us_u",D.repl="D.unemployment",
             l.patern ="ca_u|jp_u|us_u",l.repl="unemployment")

output$Canada$UECM
  
}
% \concept{ ~cpt2 }