nardl_uecm <- function(x,
                       decomp,
                       control = NULL,
                       c_q_order = NULL,
                       p_order = c(3),
                       q_order = c(4),
                       dep_var,
                       graph_save = FALSE, 
                       case = 3
){
  lagm <- function(m, nLags) #This function, lagm  is available from nardl package helper functions
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
    rhs <- cbind(dy,
                 lagy,
                 lxp,
                 lxn,
                 lagh,
                 ldy,
                 ldxp,
                 ldxn,
                 ldh,trend)
    colnames(rhs)
    rhs <- list(data = data.frame(rhs), k=k)

    y_p_order <- lagm(as.matrix(y), p_order)
    colnames(y_p_order) <- gsub(x = colnames(y_p_order), pattern = "y",replacement = dep_var)
    xp_q_order <- lagm(as.matrix(xp), q_order)
    xn_q_order <- lagm(as.matrix(xn), q_order)
    h_c_q_order <- lagm(as.matrix(h), c_q_order)
    colnames(h_c_q_order) <- gsub(x = colnames(h_c_q_order), 
                                  pattern = "x2",replacement = control)
    {
      if(p_order > 1){
        y_p_order <- y_p_order[-1,]
      }
      else{
        y_p_order <- na.omit(y_p_order)
      }
    }
    
    {
      if(c_q_order > 1){
        h_c_q_order <- h_c_q_order[-1,]
      }
      else{
        h_c_q_order <- na.omit(h_c_q_order)
      }
    }
    
    y <- data.frame(y[-1])
    colnames(y) <- dep_var
    h <- data.frame(h[-1])
    colnames(h) <- control
    nardl_data <- data.frame(cbind(y, 
                                   y_p_order,
                                   xp, xp_q_order,
                                   xn, xn_q_order,
                                   h, h_c_q_order,
                                   trend))
    nardl_case <- names(nardl_data)

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
  else {
    df <- c()
    df$y <- x[,dep_var]
    df$x1 <- x[,decomp]
    x <- model.matrix(~ x1 -1, df)
    y <- model.matrix(~ y -1, df)
    colnames(x) <- decomp
    colnames(y) <- dep_var
    lhs <- 'y'
    dy <- diff(y)
    colnames(dy) <- paste0("D.", dep_var)
    dx <- diff(x)
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
    rhs <- cbind(dy,lagy,lxp,lxn,ldy,ldxp,ldxn,trend)
    rhs <- list(data = data.frame(rhs), k=k)

        y_p_order <- lagm(as.matrix(y), p_order)
    colnames(y_p_order) <- gsub(x = colnames(y_p_order), 
                                pattern = "y",replacement = dep_var)
    xp_q_order <- lagm(as.matrix(xp), q_order)
    xn_q_order <- lagm(as.matrix(xn), q_order)
    {
      if(p_order > 1){
        y_p_order <- y_p_order[-1,]
      }
      else{
        y_p_order <- na.omit(y_p_order)
      }
    }
    
    y <- data.frame(y[-1])
    colnames(y) <- dep_var
    nardl_data <- data.frame(cbind(y, 
                                   y_p_order,
                                   xp, xp_q_order,
                                   xn, xn_q_order,
                                   trend))
    nardl_case <- names(nardl_data)
    
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
    
    dx_qlag <- lagm(as.matrix(dx), q_order)
    rhs_shr_sym <- cbind.data.frame(dy,
                                    lagy,
                                    lxp,
                                    lxn,
                                    ldy,
                                    dx_qlag,
                                    trend)
  }
  
  data <- rhs$data
  case_ <- names(data)
  {
    if(case == 1) nardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                                 paste(paste(nardl_case[-length(nardl_case)][-1], collapse= "+"),'-1')))
    if(case == 3) nardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                                 paste(nardl_case[-length(nardl_case)][-1], collapse= "+")))
    if(case == 5) nardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                                 paste(nardl_case[-1], collapse= "+")))
    if(case == 2) nardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                                 paste(nardl_case[-length(nardl_case)][-1], collapse= "+")))
    if(case == 4) nardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                                 paste(nardl_case[-1], collapse= "+")))
  }
  
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
  

  nardl_fit_case_ <- lm(nardl_fmla, data = nardl_data, na.action = na.exclude)
  fit_case_ <- lm(ecm_fmla, data = data, na.action = na.exclude)
  fit <- fit_case_
  k <- rhs$k
  
  summary_fit <- summary(fit)
  
  coeff <- summary_fit$coefficients
  nlvars <- length(coeff[,1])
  b_pos_neg <- c('_pos','_neg')
  bp <- str_subset(rownames(coeff), b_pos_neg[1])[1]
  bn <-  str_subset(rownames(coeff), b_pos_neg[2])[1]
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
      fb2 <- (-coeff[[2]]) * diag(nrow(as.matrix(fb1)))
      
    }
    else{
      fb1 <- lvars/coeff[[1]]^2
      fb2 <- (-coeff[[1]]) * diag(nrow(as.matrix(fb1)))
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
    if(is.null(control) == FALSE){
      control_lag_v <- paste0(control,'_1')
      lres <- lres[c(bp,bn,control_lag_v),]
    }
    else {
      lres <- lres[c(bp,bn),]
    }
  }

    tstat <- coef(summary(fit_case_))[paste0(dep_var, '_1'), 3]
  {
    if(!is.null(control) == T){
      if(case == 1 | case == 3 | case == 5){
        fstat_name <- c(paste0(dep_var, '_1'),bn,bp,paste0(control, '_1'))
      }
      else if(case == 2){
        fstat_name <- c("(Intercept)", paste0(dep_var, '_1'),bn,bp,paste0(control, '_1'))
      }
      else {
        fstat_name <- c("trend", paste0(dep_var, '_1'),bn,bp,paste0(control, '_1'))
      }
    }
    else {
      if(case == 1|case == 3|case == 5){
        fstat_name <- c(paste0(dep_var, '_1'),bn,bp)
      }
      else if(case == 2){
        fstat_name <- c("(Intercept)", paste0(dep_var, '_1'),bn,bp)
      }
      else {
        fstat_name <- c("trend",paste0(dep_var, '_1'),bn,bp)
      }
    }
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
  
  lr_asym_test <- cbind(asym)
  lr_asym_test <- cbind(lr_asym_test[[1]],lr_asym_test[[2]])
  colnames(lr_asym_test) <- c('statistics','p.value')
  rownames(lr_asym_test) <- decomp

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
  
  sr_asym_test <- cbind(asym)
  sr_asym_test <- cbind(sr_asym_test[[1]],sr_asym_test[[2]])
  colnames(sr_asym_test) <- c('statistics','p.value')
  rownames(sr_asym_test) <- decomp
  
  p_star <- function(x) {
    pvalue <- x[,2]
    stat <- x[,1]
    mystars <- ifelse(pvalue < .001, "*** ", ifelse(pvalue < .01, "**  ", ifelse(pvalue < .05, "*   ", ifelse(pvalue < .1, ".   ", "    "))))
    stats <- matrix(paste(stat, mystars, sep=""), nrow = nrow(x))
    message('Note: ***p < 0.001; **p < 0.01; *p < 0.05 ; (.)p < 0.1  \n')
    colnames(stats) <- 'stat'
    rownames(stats) <- rownames(x)
    return(stats)
  }
  
  jbtest <- c()
  jbtest$statistic <- round(jarque.bera.test(fit$residuals)$statistic,2)
  jbtest$p.value <- round(jarque.bera.test(fit$residuals)$p.value,2)
  jbtest <- cbind(jbtest)
  jbtest <- cbind(jbtest[[1]],jbtest[[2]])
  colnames(jbtest) <- c('JB:statistics','p.value')
  rownames(jbtest) <- decomp
  
  jbtest <- p_star(jbtest)
  
  lm_test <- c()
  lm_test$statistic <- round(bgtest(fit, type = "F", order = 4)$statistic,2)
  lm_test$p.value <- round(bgtest(fit, type = "F", order = 4)$p.value,2)
  lm_test <- cbind(lm_test)
  lm_test <- cbind(lm_test[[1]],lm_test[[2]])
  colnames(lm_test) <- c('BG:statistics','p.value')
  rownames(lm_test) <- decomp
  lm_test <- p_star(lm_test)
  
  arch <- c()
  arch$statistic <- round(ArchTest(fit$residuals, q_order)$statistic,2)
  arch$p.value <- round(ArchTest(fit$residuals, q_order)$p.value,2)
  arch <- cbind(arch)
  arch <- cbind(arch[[1]],arch[[2]])
  colnames(arch) <- c('ARCH LM:statistics','p.value')
  rownames(arch) <- decomp
  arch <- p_star(arch)
  
  diag <- cbind('lm test' = (lm_test),'arch test' = (arch), 'normality test' = (jbtest))
  colnames(diag) <- c('lm test','arch test','normality test')
  
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
  NARDL <- list('NARDL_fit' = nardl_fit_case_,
                'NARDL_ECM_fit' = fit_case_,
                'UECM' = summary_fit, 
                'cointegration'= bounds_F, 
                'Longrun_relation' = lres, 
                'longrun_asym' = lr_asym_test,
                'Shortrun_asym' = sr_asym_test, 
                'diagnostics test' = diag)
  return(NARDL)
}
