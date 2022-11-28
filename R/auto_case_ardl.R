auto_case_ardl <- function(x, dep_var, expl_var, p_order, q_order,
                           gets_pval = 0.05, order_l = 3, graph_save = FALSE){
  df <- c()
  x = na.omit(x)
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
  
  fmla <- as.formula(paste(paste0('D.', dep_var,' ~ '), 
                           paste(case_[-1], collapse= "+")))
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
  ncol(newH[[1]])
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
  
  ardl_fmla <- as.formula(paste(paste0(dep_var,' ~ '), 
                                paste(ardl_case[-1], collapse= "+")))
  
  ardl_fit_case_ <- lm(ardl_fmla, data = ardl_data, na.action = na.exclude)
  ardl_gets_fit <- gets.lm(ardl_fit_case_,
                           wald.pval=gets_pval,
                           include.1cut = T,print.searchinfo = F)
  ecm_fit_case_ <- lm(fmla, data = data, na.action = na.exclude)
  ecm_gets_fit <- gets.lm(ecm_fit_case_, wald.pval = gets_pval, keep = c(2:(1+1+length(expl_var))), include.1cut = T,print.searchinfo = F)
  ecm_gets_summ <- summary(ecm_gets_fit)
  case_6 <- ifelse(c("(Intercept)",'trend') %in% rownames(ecm_gets_summ$coefficients), 1,0)
  case = c()
  {
    if(sum(case_6) == 0){
      case = 1
    }
    if(sum(case_6) == 2){
      case = c(4,5)
    }
    
    if(sum(case_6) == 1){
      if(c("(Intercept)") %in% rownames(ecm_gets_summ$coefficients))
        case <- c(2,3)
      else if(c("trend") %in% rownames(ecm_gets_summ$coefficients))
        stop('You final model has a trend and no intercept. You should have a model with no trend and intercept; with intercept and no trend or both intercept and trend. 
      Consider adjusting the values for either the p_order, q_order or both in other to examine the various case - 1 to 5. You may consider using gets_ardl_uecm().')
    }
  }

  lrname = c()
  {
    if(length(case) < 2){
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
    }
    
    if(length(case) == 2){
      lrname <- lapply(1:length(case), function(i)
      {
        if(case[[i]] == 1| case[[i]] == 3 | case[[i]] == 5){
          lrname <- c(paste0(dep_var, '_1'),paste0(expl_var, '_1'))
        }
        else if(case[[i]] == 2){
          lrname <- c("(Intercept)", paste0(dep_var, '_1'),paste0(expl_var, '_1'))
        }
        else{
          lrname <- c("trend", paste0(dep_var, '_1'),paste0(expl_var, '_1'))
        }
      }
      )
      names(lrname) <- paste('Case', case)
    }
  }
  
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
  wld_test <- c()
  fstat <- c()
  {
    if (length(lrname) > 2){
      wld_test <- waldtest(ecm_gets_fit,test = "F", terms = lrname)
      fstat <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2], df = wld_test$Df[2])
      rownames(fstat) <- c('wald_test')
      fstat
    }
    if (length(lrname) == 2){
      for(i in 1:2){
        wld_test[[i]] <- waldtest(ecm_gets_fit,test = "F", terms = lrname[[i]])
        fstat[[i]] <- cbind(Fstat = wld_test[[i]]$F[2], Pval = wld_test[[i]]$`Pr(>F)`[2])
        rownames(fstat[[i]]) <- c('wald_test')
      }
      names(fstat) <- paste('Case', case)
    }
  }
  
  ecm_diag_gets <- round(diag <- diagnostics(ecm_gets_summ,
                                              ar.LjungB=c(order_l, 0.025), 
                                              arch.LjungB = NULL, 
                                              normality.JarqueB = NULL),3)
  
  jbtest <- c()
  jbtest$statistic <- round(jarque.bera.test(ecm_gets_fit$residuals)$statistic,2)
  jbtest$p.value <- round(jarque.bera.test(ecm_gets_fit$residuals)$p.value,2)
  jbtest <- cbind(jbtest)
  jbtest <- cbind(jbtest[[1]],jbtest[[2]])
  colnames(jbtest) <- c('statistics','p.value')
  rownames(jbtest) <- dep_var
  
  lm_test <- c()
  lm_test$statistic <- round(bgtest(ecm_gets_fit, type = "Chisq", order = order_l)$statistic,2)
  lm_test$p.value <- round(bgtest(ecm_gets_fit, type = "Chisq", order = order_l)$p.value,2)
  lm_test <- cbind(lm_test)
  lm_test <- cbind(lm_test[[1]],lm_test[[2]])
  colnames(lm_test) <- c('statistics','p.value')
  rownames(lm_test) <- dep_var
  
  arch <- c()
  arch$statistic <- round(ArchTest(ecm_gets_fit$residuals, order_l)$statistic,2)
  arch$p.value <- round(ArchTest(ecm_gets_fit$residuals, order_l)$p.value,2)
  arch <- cbind(arch)
  arch <- cbind(arch[[1]],arch[[2]])
  colnames(arch) <- c('statistics','p.value')
  rownames(arch) <- dep_var
  
  reset_test <- c()
  reset_test$statistic <- round(resettest(ecm_gets_fit, power = 2, type = 'princomp')$statistic,2)
  reset_test$p.value <- round(resettest(ecm_gets_fit, power = 2, type = 'princomp')$p.value,2)
  reset_test <- cbind(reset_test)
  reset_test <- cbind(reset_test[[1]],reset_test[[2]])
  colnames(reset_test) <- c('statistics','p.value')
  rownames(reset_test) <- dep_var
  
  diag <- rbind(lm_test, arch, jbtest, reset_test)
  rownames(diag) <- c('BG_SC_lm_test', 'LM_ARCH_test','normality_test', 'RESET_test')
  
  nobs <- nobs(ecm_gets_fit)
  stab_plot <- function(graph_save){
    oldpar <- par(no.readonly = TRUE)
    if(graph_save == TRUE){
      e <-  ecm_gets_fit$residuals
      n <-  nobs
      par(mfrow = c(1,2))
      cusum(e=e,k=k,n=n)
      cumsq(e=e,k=k,n=n)
    }
    on.exit(par(oldpar))
  }
  stab_plot(graph_save)
  
  bounds_F <- c()
  {
    if (length(case) < 2){
      {
        if(case == 1| case == 3 | case == 5){
          bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = tstat, obs=nobs,k=k)
        }
        else if(case == 2){
          bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = NULL, obs=nobs,k=k)
        }
        else{
          bounds_F <- dynamac_pkg_bounds_test(case=case,fstat=fstat[1,1],tstat = NULL, obs=nobs,k=k)
        }
      }
    }
    
    if(length(case) == 2){
      for (i in 1:2) {
        if(case[[i]] == 1| case[[i]] == 3 | case[[i]] == 5){
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case=case[[i]],fstat=fstat[[i]][1,1],tstat = tstat, obs=nobs,k=k)
        }
        else if(case[[i]] == 2){
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case=case[[i]],fstat=fstat[[i]][1,1],tstat = NULL, obs=nobs,k=k)
        }
        else{
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case=case[[i]],fstat=fstat[[i]][1,1],tstat = NULL, obs=nobs,k=k)
        }
      }
      names(bounds_F) <- paste('Case', case)
    }
  }
  
  gets_ARDL_aut <- list('Parsimonious_ARDL_fit' = ardl_gets_fit, 
                    'Parsimonious_ECM_fit' = ecm_gets_fit, 
                    'Summary_ecm_fit' = ecm_gets_summ,
                    'Parsimonious_ECM_diagnostics_test' = diag,
                    'cointegration'= bounds_F,
                    'Longrun_relation' = lres[1:dim(dh)[2],])
  
  return(gets_ARDL_aut)
}
