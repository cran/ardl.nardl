nardl_uecm_sym <- function(x,
                           assumption = c('SRSR'),
                           decomp,
                           control = NULL,
                           c_q_order = NULL,
                           p_order = c(3),
                           q_order = c(4),
                           dep_var,
                           graph_save = FALSE, 
                           case = 3){
  lagm <- function(m, nLags)
  {
    if(!is.matrix(m))
      stop("Trying to lag something that's not a matrix")
    
    d <- dim(m)
    if(is.null(colnames(m)))
      colnames(m) <- as.character(seq_len(d[2]))
    if(nLags > d[1])
      stop(sprintf("You try to create %d lags but there's only %d rows in m.",
                   nLags, d[1]))
    
    lagM <- matrix(NA,nrow=d[1], ncol = d[2]*nLags)
    
    for(i in seq_len(nLags)){
      cid <- seq(1:d[2]) + d[2]*(i-1)
      
      lagM[(i+1):d[1],cid] <- m[1:(d[1]-i),]
    }
    
    cnames <- outer(colnames(m),seq_len(nLags), FUN = paste, sep = "_")
    
    colnames(lagM) <- c(cnames)
    
    return(lagM)
    
  }
  if(assumption %in% c( 'LRSR','SRSR') == FALSE)
    warning('assumption should take either LRSR or SRSR.')
  {
    if(!is.null(control) == T){
      df <- c()
      df$y <- x[,dep_var]
      df$x1 <- x[,decomp]
      df$x2 <- x[,control]
      x <- model.matrix(~ x1 -1, df)
      y <- model.matrix(~ y -1, df)
      colnames(x) <- decomp
      dy <- diff(y)
      colnames(dy) <- paste("D", dep_var, sep = ".")
      dx <- diff(x)
      colnames(dx) <- paste("D", colnames(x), sep = ".")
      h <- model.matrix(~ x2 -1, df)
      dh <- diff(h)
      colnames(dh) <- paste0("D.", control)
      n <- nrow(dx)
      cl <- ncol(dx)
      pos <- dx[, 1:cl] >= 0
      dxp <- as.matrix(as.numeric(pos) * dx[, 1:cl])
      colnames(dxp) <- paste(colnames(dx), "pos", sep = "_")
      dxn <- as.matrix((1 - as.numeric(pos)) * dx[, 1:cl])
      colnames(dxn) <- paste(colnames(dx), "neg", sep = "_")
      xp <- apply(dxp, 2, cumsum)
      colnames(xp) <- paste(colnames(x), "pos", sep = "_")
      xn <- apply(dxn, 2, cumsum)
      colnames(xn) <- paste(colnames(x), "neg", sep = "_")
      lagy <- lagm(as.matrix(y), 1)
      colnames(lagy) <- paste0(dep_var,'_1')
      lagh <- lagm(as.matrix(h), 1)
      colnames(lagh) <- paste0(control,'_1')
      lxp <- lagm(as.matrix(xp), 1)
      lxn <- lagm(as.matrix(xn), 1)
      ldy <- lagm(dy, p_order)
      ldh <- lagm(as.matrix(dh), c_q_order)
      ldxp <- lagm(as.matrix(dxp), q_order)
      ldxn <- lagm(as.matrix(dxn), q_order)
      lagy <- na.omit(lagy)
      lagh <- na.omit(lagh)
      trend <- seq_along(dy)
      k <- ncol(x) + ncol(h)
      if(assumption == 'LRSR'){
        lagx <- lagm(as.matrix(x), 1)
        lagx1 <- data.frame(lagx[-1])
        colnames(lagx1) <- colnames(lagx)
        rhs_lng_sym <- cbind(dy,
                             lagy,
                             lagx1,
                             lagh,
                             ldy,
                             ldxp,
                             ldxn,
                             ldh,
                             trend)
        case_ <- names(rhs_lng_sym)
      }
      else if(assumption == 'SRSR'){
        dx_qlag <- lagm(as.matrix(dx), q_order)
        rhs_shr_sym <- cbind.data.frame(dy,
                                        lagy,
                                        lxp,
                                        lxn,
                                        lagh,
                                        ldy,
                                        dx_qlag,
                                        ldh,
                                        trend)
      }
      case_ <- names(rhs_shr_sym)
    }
  }
  
  {
    if(!is.null(control) == F){
      df <- c()
      df$y <- x[,dep_var]
      df$x1 <- x[,decomp]
      (x <- model.matrix(~ x1 -1, df))
      (y <- model.matrix(~ y -1, df))
      colnames(x) <- decomp
      colnames(y) <- dep_var
      lhs <- 'y'
      (dy <- diff(y))
      colnames(dy) <- paste0("D.", dep_var)
      (dx <- diff(x))
      colnames(dx) <- paste0("D.", colnames(x))
      n <- nrow(dx)
      cl <- ncol(dx)
      pos <- dx[, 1:cl] >= 0
      dxp <- as.matrix(as.numeric(pos) * dx[, 1:cl])
      colnames(dxp) <- paste(colnames(dx), "pos", sep = "_")
      dxn <- as.matrix((1 - as.numeric(pos)) * dx[, 1:cl])
      colnames(dxn) <- paste(colnames(dx), "neg", sep = "_")
      xp <- apply(dxp, 2, cumsum)
      colnames(xp) <- paste(colnames(x), "pos", sep = "_")
      xn <- apply(dxn, 2, cumsum)
      colnames(xn) <- paste(colnames(x), "neg", sep = "_")
      lagy <- lagm(as.matrix(y), 1)
      lxp <- lagm(as.matrix(xp), 1)
      lxn <- lagm(as.matrix(xn), 1)
      k <- ncol(x)
      ldy <- lagm(dy, p_order)
      ldxp <- lagm(as.matrix(dxp), q_order)
      ldxn <- lagm(as.matrix(dxn), q_order)
      lagy <- na.omit(lagy)
      trend <- seq_along(dy)
      if(assumption == 'LRSR'){
        lagx <- lagm(as.matrix(x), 1)
        lagx1 <- data.frame(lagx[-1])
        colnames(lagx1) <- colnames(lagx)
        rhs_lng_sym <- cbind(dy,
                             lagy,
                             lagx1,
                             ldy,
                             ldxp,
                             ldxn,
                             trend)
        case_ <- names(rhs_lng_sym)
      }
      if(assumption == 'SRSR'){
        dx_qlag <- lagm(as.matrix(dx), q_order)
        rhs_shr_sym <- cbind.data.frame(dy,
                                        lagy,
                                        lxp,
                                        lxn,
                                        ldy,
                                        dx_qlag,
                                        trend)
        case_ <- names(rhs_shr_sym)
      }
    }
  }

  {
    if(assumption == 'LRSR'){
      case_ <-names(rhs_lng_sym)
      data <- rhs_lng_sym
      {
        if(case == 1) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(paste(case_[-length(case_)][-1], collapse= "+"),'-1')))
        if(case == 3) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-length(case_)][-1], collapse= "+")))
        if(case == 5) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-1], collapse= "+")))
        if(case == 2) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-length(case_)][-1], collapse= "+")))
        if(case == 4) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-1], collapse= "+")))
      }
    }
    if(assumption == 'SRSR'){
      case_ <- names(rhs_shr_sym)
      data <- rhs_shr_sym
      {
        if(case == 1) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(paste(case_[-length(case_)][-1], collapse= "+"),'-1')))
        if(case == 3) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-length(case_)][-1], collapse= "+")))
        if(case == 5) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-1], collapse= "+")))
        if(case == 2) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-length(case_)][-1], collapse= "+")))
        if(case == 4) ecm_fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                                                   paste(case_[-1], collapse= "+")))
      }
    }
    ecm_fmla
  }
  fit_case_ <- lm(ecm_fmla, data = data, na.action = na.exclude)
  fit <- fit_case_
  summary_fit <- summary(fit)
  
  coeff <- summary_fit$coefficients
  nlvars <- length(coeff[,1])
  b_pos_neg <- c('_pos','_neg')
  
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
  vc <- vcov(fit)
  {
    if("(Intercept)" %in% rownames(coeff) == TRUE){
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

    {
    if(assumption == 'LRSR'){
      if(!is.null(control) == T){
        control_lag_v <- str_subset(rownames(coeff), paste0(control,'_1'))[1]
        decomp_lag_v <- str_subset(rownames(coeff), paste0(decomp,'_1'))
        lres <- lres[c(decomp_lag_v,control_lag_v),]
      }
      else {
        decomp_lag_v <- str_subset(rownames(coeff), paste0(decomp,'_1'))
        lres <- data.frame(lres)
        lres <- lres[c(decomp_lag_v),]
      }
    }
    
    if(assumption == 'SRSR'){
      if(!is.null(control) == T){
        bp <- str_subset(rownames(coeff), b_pos_neg[1])
        bn <-  str_subset(rownames(coeff), b_pos_neg[2])
        control_lag_v <- str_subset(rownames(coeff), paste0(control,'_1'))[1]
        lres <- lres[c(bp,bn,control_lag_v),]
      }
      else {
        bp <- str_subset(rownames(coeff), b_pos_neg[1])
        bn <-  str_subset(rownames(coeff), b_pos_neg[2])
        lres <- lres[c(bp,bn),]
      }
    }
    lres
  }
  
  tstat <- coef(summary(fit_case_))[paste0(dep_var, '_1'), 3]
  fstat_name <- rownames(lres)
  {
    if(case == 1 | case == 3 | case == 5){
      fstat_name <- c(paste0(dep_var, '_1'),fstat_name)
    }
    else if(case == 2){
      fstat_name <- c("(Intercept)", paste0(dep_var, '_1'),fstat_name)
    }
    else {
      fstat_name <- c("trend", paste0(dep_var, '_1'),fstat_name)
    }
    fstat_name
  }
  
  fstat_name <- paste0(fstat_name, "=0")
  coin <- linearHypothesis(fit, fstat_name, test = "F")
  fstat <- coin$F[2]
  {
    if("(Intercept)" %in% rownames(vc) == TRUE){
      tcoin <- fit_case_$coefficients[2]
    }
    else{
      tcoin <- fit_case_$coefficients[1]
    }
  }
  
  {
    if(assumption == 'LRSR'){
      {
        decomp_d = paste0('D.', decomp)
        sr_cof_nm <- str_subset(rownames(coeff), decomp_d)
        coeff_sr <- coeff[sr_cof_nm,]
        
        d_pos <-  sr_cof_nm[1:(length(sr_cof_nm)/2)]
        d_neg <-  sr_cof_nm[(length(d_pos)+1):length(sr_cof_nm)]
        
        dddneg <-  paste(d_neg,'+')
        dddpos <-  paste(d_pos,'+')
        
        plus_ <-  dddpos[[length(dddpos)]]
        plus_ <-  gsub(x = plus_, pattern = "[+]| ",replacement = "")
        plus_ <-  paste(plus_,'=')
        
        dddpos <- str_flatten(c(dddpos[-length(dddpos)],plus_))
        
        neg_ <-  dddneg[[length(dddneg)]]
        neg_ <-  gsub(x = neg_, pattern = "[+]| ",replacement = "")
        dddneg <- str_flatten(c(dddneg[-length(dddneg)],neg_))
        
        pos_neg <- str_c(dddpos, dddneg)
        
        asym <- c()
        asym$F <-  linearHypothesis(fit, pos_neg, test = "F",white.adjust = F )$F[2]
        
        asym$`Pr(>F)` <- linearHypothesis(fit, pos_neg, test = "F",white.adjust = F )$`Pr(>F)`[2]
        
        asym_test <- cbind(asym)
        asym_test <- cbind(asym_test[[1]],asym_test[[2]])
        colnames(asym_test) <- c('statistics','p.value')
        rownames(asym_test) <- paste0('Shortrun asymmetric:', decomp)    }
    }
    if(assumption == 'SRSR'){
      {
        rcap = matrix(c(0), 1, 2)
        rcap[1, 1] = 1
        rcap[1, 2] = -1
        rcap
        rsml = matrix(c(0), 1, 1)
        coof_lres <- coof[c(bn,bp)]
        vcc1 <- vcc[c(names(coof_lres[1]), names(coof_lres[2])), 
                     c(names(coof_lres[1]), names(coof_lres[2]))]
        
        lbeta <- coeff[c(names(coof_lres[1]), names(coof_lres[2])), 1]
        ll2 <- c(paste(names(coof_lres[1]), "=", names(coof_lres[2])))
        
        asym <- c()
        asym$F <- linearHypothesis(fit, ll2, test = "F",white.adjust = F )$F[2]
        
        asym$`Pr(>F)` <- linearHypothesis(fit, ll2, test = "F",white.adjust = F )$`Pr(>F)`[2]
        
        asym_test <- cbind(asym)
        asym_test <- cbind(asym_test[[1]],asym_test[[2]])
        colnames(asym_test) <- c('statistics','p.value')
        rownames(asym_test) <- paste0('Longrun asymmetric:', decomp)
        
      }
    }
    asym_test
  }
  
  jbtest <- c()
  jbtest$statistic <- round(jarque.bera.test(fit$residuals)$statistic,3)
  jbtest$p.value <- round(jarque.bera.test(fit$residuals)$p.value,3)
  jbtest <- cbind(jbtest)
  jbtest <- cbind(jbtest[[1]],jbtest[[2]])
  colnames(jbtest) <- c('JB:statistics','p.value')
  rownames(jbtest) <- decomp
  
  lm_test <- c()
  lm_test$statistic <- round(bgtest(fit, type = "F", order = 4)$statistic,3)
  lm_test$p.value <- round(bgtest(fit, type = "F", order = 4)$p.value,3)# has to be > 0.05 
  lm_test <- cbind(lm_test)
  lm_test <- cbind(lm_test[[1]],lm_test[[2]])
  colnames(lm_test) <- c('BG:statistics','p.value')
  rownames(lm_test) <- decomp
  
  arch <- c()
  arch$statistic <- round(ArchTest(fit$residuals, q_order)$statistic,2)
  arch$p.value <- round(ArchTest(fit$residuals, q_order)$p.value,2)
  arch <- cbind(arch)
  arch <- cbind(arch[[1]],arch[[2]])
  colnames(arch) <- c('ARCH LM:statistics','p.value')
  rownames(arch) <- decomp
  
  diag <- cbind('lm test' = (lm_test),'arch test' = (arch), 'normality test' = (jbtest))
  colnames(diag) <- c('lm test','pvalue', 'arch test','pvalue','normality test','pvalue')
  diag
  nobs <- nobs(fit)
  e <-  fit$residuals
  stab_plot <- function(graph_save){
    oldpar <- par(no.readonly = TRUE)
    if(graph_save == TRUE){
      e <-  fit$residuals
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
  
  NARDL <- list(
    'NARDL_ECM_fit' = fit_case_,
    'UECM' = summary_fit, 
    'cointegration'= bounds_F, 
    'Longrun_relation' = lres, 
    'asymmetric test' = asym_test,
    'diagnostics test' = diag)
  return(NARDL)
}
