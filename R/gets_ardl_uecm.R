gets_ardl_uecm <- function(x, dep_var, expl_var, p_order = c(2), 
                    q_order = c(3), gets_pval = 0.1, case = 3,
                    F_HC = FALSE, order_l = 5, graph_save = FALSE){
df <- c()
x <- na.omit(x)
df$y <- x[,dep_var]
y <- model.matrix(~ y -1, df)
lhs <- 'y'
dy <- diff(y)
colnames(dy) <- paste("D", dep_var, sep = ".")
{
  if(length(expl_var > 1)){
    expl2 <- paste(paste(expl_var,'+'))[-length(expl_var)]
    expl3 <- paste(paste(expl_var,'+'))[length(expl_var)]
    expl3 <- gsub(x = expl3, pattern = "[+]",replacement = "-1")
    ff <- str_flatten(c(expl2,expl3))
    ff <- formula(paste('~',ff))
  }
  else{
    expl <- paste(paste(expl_var,'+'))
    expl1 <- gsub(x = expl, pattern = "[+]",replacement = "-1")
    ff <- str_flatten(c(expl1))
    ff <- formula(paste('~',ff))
  }
}
h <- model.matrix(ff, x)
dh <- diff(h)
colnames(dh) <- paste("D", colnames(h), sep = ".")
lagy <- lagm(y, 1)
colnames(lagy) <- paste0(dep_var,'_1')
lagh <- lagm((h), 1)
k <- ncol(h)
ldy <- lagm(dy, p_order)
ldh <- lapply(1:dim(dh)[2], function(i) 
  lagm(as.matrix(dh[,i]), q_order[i]))
colnames_ldh <- lapply(1:dim(dh)[2], function(i)
  gsub(x = colnames(ldh[[i]]),
       pattern = c('_'),
       replacement = paste0(colnames(dh)[i],'_')))
for (i in 1:dim(dh)[2]) {
  colnames(ldh[[i]]) <- colnames_ldh[[i]]
}
lagy <- na.omit(lagy)
lagh <- na.omit(lagh)
trend <- seq_along(dy)
NCOL <-lapply(1:dim(dh)[2], function(i)
  ncol(ldh[[i]]))
ldh <- data.frame(matrix(unlist(ldh),
                         ncol = sum(unlist(NCOL)), byrow = F))
colnames(ldh) <- unlist(colnames_ldh)
rhnams <- c("Const", 
            colnames(dy),
            colnames(lagy), 
            colnames(lagh),
            colnames(ldy), 
            colnames(dh),
            colnames(ldh),'trend')
data_uecm <- cbind(dy,
                   lagy,
                   lagh, 
                   ldy,
                   dh,
                   ldh,
                   trend)
data_uecm <- list(data = data.frame(data_uecm), k=k)
data <- data_uecm$data
case_ <- names(data)

{
  if(case == 1) fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                         paste(paste(case_[-length(case_)][-1], collapse= "+"),'-1')))
  if(case == 3) fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                         paste(case_[-length(case_)][-1], collapse= "+")))
  if(case == 5) fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                         paste(case_[-1], collapse= "+")))
  if(case == 2) fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                         paste(case_[-length(case_)][-1], collapse= "+")))
  if(case == 4) fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                         paste(case_[-1], collapse= "+")))
}
y_p_order <- lagm((y), p_order)
colnames(y_p_order) <- gsub(x = colnames(y_p_order), pattern = "y",
                            replacement = dep_var)
newH <- lapply(1:dim(h)[2], function(i)
  lagm(as.matrix(h[,i]), q_order[[i]]))
colnames_newH <- lapply(1:dim(h)[2], function(i)
  gsub(x = colnames(newH[[i]]),
       pattern = c('_'),
       replacement = paste0(colnames(h)[i],'_')))
for (i in 1:dim(h)[2]) {
  colnames(newH[[i]]) <- colnames_newH[[i]]
}

NCOL <-lapply(1:dim(h)[2], function(i) 
  ncol(newH[[i]]))
newH <- data.frame(matrix(unlist(newH),
                          ncol = sum(unlist(NCOL)), byrow = F))
colnames(newH) <- unlist(colnames_newH)
trend <- seq(y_p_order[,1])
colnames(y) <- dep_var
ardl_data <- cbind(y,
                   y_p_order,
                   h,
                   newH, 
                   trend)
ardl_case <- names(ardl_data)
{
  if(case == 1) ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                              paste(paste(ardl_case[-length(ardl_case)][-1], collapse= "+"),'-1')))
  if(case == 3) ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                              paste(ardl_case[-length(ardl_case)][-1], collapse= "+")))
  if(case == 5) ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                              paste(ardl_case[-1], collapse= "+")))
  if(case == 2) ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                              paste(ardl_case[-length(ardl_case)][-1], collapse= "+")))
  if(case == 4) ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                              paste(ardl_case[-1], collapse= "+")))
}

ardl_fit_case_ <- lm(ardl_fmla, data = ardl_data, na.action = na.exclude)

ardl_gets_fit <- gets.lm(ardl_fit_case_,
                         wald.pval=gets_pval,
                         include.1cut = T,print.searchinfo = F)

ecm_fit_case_ <- lm(fmla, data = data, na.action = na.exclude)

{
  if("(Intercept)" %in% names(ecm_fit_case_$coefficients) == TRUE){
    keep_val <- c(1:(1+1+length(expl_var)))
  }
  else{
    keep_val <- c(1:(1+length(expl_var)))
  }
}

{
  if("trend" %in% names(ecm_fit_case_$coefficients) == TRUE){
    keep_val <- c(keep_val,length(ecm_fit_case_$coefficients))
  }
  else{
    keep_val <- keep_val
  }
}

ecm_gets_fit <- gets.lm(ecm_fit_case_,
                        wald.pval=gets_pval,
                        keep = keep_val,
                        include.1cut = T,print.searchinfo = F)
ecm_gets_summ <- summary(ecm_gets_fit)

coeff <- ecm_gets_summ$coefficients
nlvars <- length(coeff[,1])
{
  if("(Intercept)" %in% rownames(coeff) == TRUE){
    lvars <- coeff[3:nlvars, 1]
    coof <- -lvars/coeff[[2]]
  }
  else{
    lvars <- coeff[2:nlvars, 1]
    coof <- -lvars/coeff[[1]]
  }
}

seldata <- data.matrix(coeff)
cof <- matrix(coof, length(lvars), 1)

{
  if("(Intercept)" %in% rownames(coeff) == TRUE){
    fb1 <- lvars/coeff[[2]]^2
    fb2 <- (-1/coeff[[2]]) * diag(nrow(as.matrix(fb1)))
    
  }
  else{
    fb1 <- lvars/coeff[[1]]^2
    fb2 <- (-1/coeff[[1]]) * diag(nrow(as.matrix(fb1)))
  }
}

fb <- cbind(as.matrix(fb1), fb2)
vc <- vcov(ecm_gets_fit)

{
  if("(Intercept)" %in% rownames(vc) == TRUE){
    vcc <-vc[2:nrow(vc), 2:ncol(vc)]
  }
  else{
    vcc <-vc[1:nrow(vc), 1:ncol(vc)]
  }
}

lrse <- sqrt(diag(fb %*% vcc %*% t(fb)))
lrt <- coof/lrse
X <- sum(ifelse(names(ecm_gets_fit$model) %in% names(coef(ecm_gets_fit)), 1, 0))
Y <- dim(ecm_gets_fit$model)[1]
rdf <- Y - X
lrpv <- 2*pt(abs(coof/lrse),df = rdf,lower.tail=FALSE)
lres <- cbind(coof,lrse, lrt,lrpv)
colnames(lres) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

tstat <- coef(summary(ecm_gets_fit))[paste0(dep_var, '_1'), 3]
lrname <- c()
{
  if(case == 1| case == 3 | case == 5){
    lrname <- c(paste0(dep_var, '_1'),paste0(expl_var, '_1'))
  }
  else if(case == 2){
    lrname <- c("(Intercept)", paste0(dep_var, '_1'),paste0(expl_var, '_1'))
  }
  else{
    lrname <- c("trend", paste0(dep_var, '_1'),paste0(expl_var, '_1'))
  }
}

wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = lrname, verbose = F,
                              rhs = rep(0,length(lrname)), test = "F")
wld_test_ <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = lrname, verbose = F,
                               rhs = rep(0,length(lrname)), test = "F",vcov = vcovHC)
lmwld <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2], df = wld_test$Df[2])
lmwld_ <- cbind(Fstat = wld_test_$F[2], Pval = wld_test_$`Pr(>F)`[2], df = wld_test_$Df[2])

fstat <- rbind(lmwld,lmwld_)
rownames(fstat) <- c('wald_test','robust_wald_test')
ecm_diag_gets <- round(diag <- diagnostics(ecm_gets_summ,
                                            ar.LjungB=c(floor(max(q_order)/2), 0.025), 
                                            arch.LjungB = NULL, 
                                            normality.JarqueB = NULL),3)


jbtest <- c()
jbtest$statistic <- round(jarque.bera.test(ecm_gets_fit$residuals)$statistic,4)
jbtest$p.value <- round(jarque.bera.test(ecm_gets_fit$residuals)$p.value,4)
jbtest <- cbind(jbtest)
jbtest <- cbind(jbtest[[1]],jbtest[[2]])
colnames(jbtest) <- c('statistics','p.value')
rownames(jbtest) <- dep_var

lm_test <- c()
lm_test$statistic <- round(bgtest(ecm_gets_fit, type = "Chisq", order = order_l)$statistic,4)
lm_test$p.value <- round(bgtest(ecm_gets_fit, type = "Chisq", order = order_l)$p.value,4)
lm_test <- cbind(lm_test)
lm_test <- cbind(lm_test[[1]],lm_test[[2]])
colnames(lm_test) <- c('statistics','p.value')
rownames(lm_test) <- dep_var

arch <- c()
arch$statistic <- round(ArchTest(ecm_gets_fit$residuals, order_l)$statistic,4)
arch$p.value <- round(ArchTest(ecm_gets_fit$residuals, order_l)$p.value,4)
arch <- cbind(arch)
arch <- cbind(arch[[1]],arch[[2]])
colnames(arch) <- c('statistics','p.value')
rownames(arch) <- dep_var

reset_test <- c()
reset_test$statistic <- round(resettest(ecm_gets_fit, power = 2, type = 'princomp')$statistic,4)
reset_test$p.value <- round(resettest(ecm_gets_fit, power = 2, type = 'princomp')$p.value,4)
reset_test <- cbind(reset_test)
reset_test <- cbind(reset_test[[1]],reset_test[[2]])
colnames(reset_test) <- c('statistics','p.value')
rownames(reset_test) <- dep_var

diag <- rbind(lm_test, arch, jbtest, reset_test)
rownames(diag) <- c('BG_SC_lm_test', 'LM_ARCH_test','normality_test', 'RESET_test')

nobs <- nobs(ecm_gets_fit)
stab_plot <- function(graph_save){
  if(graph_save == TRUE){
    oldpar <- par(no.readonly = TRUE)
    e <-  ecm_gets_fit$residuals
    n <-  nobs
    par(mfrow = c(1,2))
    cusum(e=e,k=k,n=n)
    cumsq(e=e,k=k,n=n)
    on.exit(par(oldpar))
  }
}
stab_plot(graph_save)

{
  if(F_HC == T)
    if(case == 1| case == 3 | case == 5){
      if(diag[2,2] < 0.05)
        bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[2,1],tstat = tstat, obs=nobs,k=k)
      else if (diag[2,2] > 0.05)
        bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = tstat, obs=nobs,k=k)
    }
  else if(case == 2){
    if(diag[2,2] < 0.05)
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[2,1],tstat = NULL, obs=nobs,k=k)
    else if (diag[2,2] > 0.05)
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = NULL, obs=nobs,k=k)
  }
  else{
    if(diag[2,2] < 0.05)
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[2,1],tstat = NULL, obs=nobs,k=k)
    else if (diag[2,2] > 0.05)
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = NULL, obs=nobs,k=k)
  }
  else{
    if (case == 1 | case == 3 | case == 5) {
      bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                          fstat = fstat[1,1], tstat = tstat, obs = nobs, k = k)
    }
    else if (case == 2) {
      bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                          fstat = fstat[1,1], tstat = NULL, obs = nobs, k = k)
    }
    else {
      bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                          fstat = fstat[1,1], tstat = NULL, obs = nobs, k = k)
    }
  }
}

diag <- rbind(diag,ecm_diag_gets[,-2])
rownames(diag) <- c(rownames(diag)[-5],rownames(ecm_diag_gets))

gets_ARDL <- list('Parsimonious_ARDL_fit' = ardl_gets_fit, 
                  'Parsimonious_ECM_fit' = ecm_gets_fit, 
                  'Summary_ecm_fit' = ecm_gets_summ,
                  'Parsimonious_ECM_diagnostics_test' = diag,
                  'cointegration'= bounds_F, 
                  'Longrun_relation' = lres[1:dim(dh)[2],])
return(gets_ARDL)
}
