ardl_uecm <- function(x,
               p_order = c(2), 
               q_order,
               dep_var,
               order_l = 4,
               graph_save = FALSE,
               expl_var,
               case = 3){
  
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
  lagy <- lagm(as.matrix(y), 1)
  colnames(lagy) <- paste0(dep_var,'_1')
  lagh <- lagm(as.matrix(h), 1)
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
  rhnams
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
  
  y_p_order <- lagm(as.matrix(y), p_order)
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
  fit_case_ <- lm(fmla, data = data, na.action = na.exclude)
  k <- data_uecm$k
  summary_fit <- summary(fit_case_)
  
  coeff <- summary_fit$coefficients
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
  vc <- vcov(fit_case_)
  
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
  X <- sum(ifelse(names(fit_case_$model) %in% names(coef(fit_case_)), 1, 0))
  Y <- dim(fit_case_$model)[1]
  rdf <- Y - X
  lrpv <- 2*pt(abs(coof/lrse),df = rdf,lower.tail=FALSE)
  lres <- cbind(coof,lrse, lrt,lrpv)
  
  colnames(lres) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  tstat <- coef(summary(fit_case_))[paste0(dep_var, '_1'), 3]
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
  l <- paste0(lrname, "=0")
  coin <- linearHypothesis(fit_case_, l, test = "F")
  
  fstat <- coin$F[2]
  {
    if("(Intercept)" %in% rownames(vc) == TRUE){
      tcoin <- fit_case_$coefficients[2]
    }
    else{
      tcoin <- fit_case_$coefficients[1]
    }
  }
  
  jbtest <- c()
  jbtest$statistic <- round(jarque.bera.test(fit_case_$residuals)$statistic,4)
  jbtest$p.value <- round(jarque.bera.test(fit_case_$residuals)$p.value,4)
  jbtest <- cbind(jbtest)
  jbtest <- cbind(jbtest[[1]],jbtest[[2]])
  colnames(jbtest) <- c('statistics','p.value')
  rownames(jbtest) <- dep_var
  
  lm_test <- c()
  lm_test$statistic <- round(bgtest(fit_case_, type = "Chisq", order = order_l)$statistic,4)
  lm_test$p.value <- round(bgtest(fit_case_, type = "Chisq", order = order_l)$p.value,4)
  lm_test <- cbind(lm_test)
  lm_test <- cbind(lm_test[[1]],lm_test[[2]])
  colnames(lm_test) <- c('statistics','p.value')
  rownames(lm_test) <- dep_var
  
  arch <- c()
  arch$statistic <- round(ArchTest(fit_case_$residuals, order_l)$statistic,4)
  arch$p.value <- round(ArchTest(fit_case_$residuals, order_l)$p.value,4)
  arch <- cbind(arch)
  arch <- cbind(arch[[1]],arch[[2]])
  colnames(arch) <- c('statistics','p.value')
  rownames(arch) <- dep_var
  
  reset_test <- c()
  reset_test$statistic <- round(resettest(fit_case_, power = 2, type = 'princomp')$statistic,4)
  reset_test$p.value <- round(resettest(fit_case_, power = 2, type = 'princomp')$p.value,4)
  reset_test <- cbind(reset_test)
  reset_test <- cbind(reset_test[[1]],reset_test[[2]])
  colnames(reset_test) <- c('statistics','p.value')
  rownames(reset_test) <- dep_var
  
  diag <- rbind(lm_test, arch, jbtest, reset_test)
  rownames(diag) <- c('BG_SC_lm_test', 'LM_ARCH_test','normality_test', 'RESET_test')
  
  nobs <- nobs(fit_case_)
  stab_plot <- function(graph_save){
    oldpar <- par(no.readonly = TRUE)
    if(graph_save == TRUE){
      e <-  fit_case_$residuals
      n <-  nobs
      par(mfrow = c(1,2))
      cusum(e=e,k=k,n=n)
      cumsq(e=e,k=k,n=n)
    }
    on.exit(par(oldpar))
  }
  stab_plot(graph_save)
  
  {
    if(case == 1| case == 3 | case == 5){
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat,tstat = tstat, obs=nobs,k=k)
    }
    else if(case == 2){
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat,tstat = NULL, obs=nobs,k=k)
    }
    else{
      bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat,tstat = NULL, obs=nobs,k=k)
    }
  }
  
  ARDL <- list('ARDL_fit' = ardl_fit_case_, 
               'ARDL_ECM_fit' = fit_case_, 
               'UECM' = summary_fit, 
               'cointegration'= bounds_F, 
               'Longrun_relation' = lres[1:dim(dh)[2],], 
               'diagnostics test' = diag)
  return(ARDL)
}
