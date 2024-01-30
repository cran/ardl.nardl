nardl_mdv <- function(x,
                      dep_var,
                      decomp1, decomp2,
                      thresh1 = Inf,
                      thresh2 = Inf,
                      gets = TRUE,
                      gets_pval = 0.1,
                      case = NULL,
                      conservative = FALSE,
                      p_order = c(3), 
                      q_order1 = c(5),
                      q_order2 = c(5),
                      order_l = 4,
                      graph_save = FALSE){
  {
    x <- na.omit(x)
    df <- c()
    df$y <- x[, dep_var]
    df$x1 <- x[, decomp1]
    df$x2 <- x[, decomp2]
    x1 <- model.matrix(~x1 - 1, df)
    x2 <- model.matrix(~x2 - 1, df)
    y <- model.matrix(~y - 1, df)
    colnames(x1) <- decomp1
    colnames(x2) <- decomp2
    dy <- diff(y)
    colnames(dy) <- paste("D", dep_var, sep = ".")
    dx1 <- diff(x1)
    dx2 <- diff(x2)
    colnames(dx1) <- paste("D", colnames(x1), sep = ".")
    colnames(dx2) <- paste("D", colnames(x2), sep = ".")
    
    n1 <- nrow(dx1)
    n2 <- nrow(dx2)
    cl <- ncol(dx1)
    
    pos1 <- dx1[, 1:cl] >= 0
    pos2 <- dx2[, 1:cl] >= 0
    
    dxp1 <- as.matrix(as.numeric(pos1) * dx1[, 1:cl])
    dxp2 <- as.matrix(as.numeric(pos2) * dx2[, 1:cl])
    
    colnames(dxp1) <- paste(colnames(dx1), "pos", sep = "_")
    colnames(dxp2) <- paste(colnames(dx2), "pos", sep = "_")
    
    dxn1 <- as.matrix((1 - as.numeric(pos1)) * dx1[, 1:cl])
    dxn2 <- as.matrix((1 - as.numeric(pos2)) * dx2[, 1:cl])
    
    colnames(dxn1) <- paste(colnames(dx1), "neg", sep = "_")
    colnames(dxn2) <- paste(colnames(dx2), "neg", sep = "_")
    
    message(ifelse(ceiling((sum(pos1)/(n1-1))*100) == 50, 
                   'positive and negative change in decomp1 are equal.', 
                   paste('Percentage of positive change in decomp1 is',
                         ceiling((sum(pos1)/(n1-1))*100), 'percent while negative change is',
                         100 - ceiling((sum(pos1)/(n1-1))*100))))
    
    message(ifelse(ceiling((sum(pos2)/(n2-1))*100) == 50, 
                   'positive and negative changes in decomp2 are equal.', 
                   paste('Percentage of positive changes in decomp2 is',
                         ceiling((sum(pos2)/(n2-1))*100), 'percent while negative change is',
                         100 - ceiling((sum(pos2)/(n2-1))*100))))
    
    cumsum_reset <- function(threshold) {
      function(x) {
        accumulate(x, ~if_else(.x>=threshold, .y, .x+.y))
      }  
    }
    
    nacumsum1 <- function(x, ...) c(cumsum_reset(threshold = thresh1[[1]])(x, ...))
    nacumsum2 <- function(x, ...) c(cumsum_reset(threshold = thresh2[[1]])(x, ...))
    
    if (thresh1 == 'mean'){
      (thresh1 = mean(dx1))
    }
    
    if (thresh2 == 'mean'){
      (thresh2 = mean(dx2))
    }
    
    xp1 <- apply(dxp1, 2, nacumsum1)
    xp2 <- apply(dxp2, 2, nacumsum2)
    
    colnames(xp1) <- paste(colnames(x1), "pos", sep = "_")
    colnames(xp2) <- paste(colnames(x2), "pos", sep = "_")
    
    xn1 <- apply(dxn1, 2, nacumsum1)
    xn2 <- apply(dxn2, 2, nacumsum2)
    
    colnames(xn1) <- paste(colnames(x1), "neg", sep = "_")
    colnames(xn2) <- paste(colnames(x2), "neg", sep = "_")
    
    lagy <- lagm(as.matrix(y), 1)
    colnames(lagy) <- paste0(dep_var, "_1")
    
    lxp1 <- lagm(as.matrix(xp1), 1)
    lxp2 <- lagm(as.matrix(xp2), 1)
    
    lxn1 <- lagm(as.matrix(xn1), 1)
    lxn2 <- lagm(as.matrix(xn2), 1)
    
    ldy <- lagm(dy, p_order)
    
    ldxp1 <- lagm(as.matrix(dxp1), q_order1)
    ldxp2 <- lagm(as.matrix(dxp2), q_order2)
    
    ldxn1 <- lagm(as.matrix(dxn1), q_order1)
    ldxn2 <- lagm(as.matrix(dxn2), q_order2)
    
    lagy <- na.omit(lagy)
    trend <- seq_along(dy)
    
    {
      if(conservative == TRUE){
        (k <- ncol(x1) +ncol(x2))
      }
      else{
        (k <- ncol(x1) + ncol(x2) + 2)
      }
    }
    
    rhs <- cbind(dy, lagy, 
                 lxp1, lxp2, 
                 lxn1, lxn2, 
                 ldy, 
                 ldxp1, ldxp2, 
                 ldxn1, ldxn2,
                 trend)
    
    rhs <- cbind(dy, lagy, 
                 lxp1, lxn1,
                 lxp2, lxn2, 
                 ldy, 
                 ldxp1, ldxn1,
                 ldxp2, ldxn2,
                 trend)
    colnames(rhs)
    rhs <- list(data = data.frame(rhs), k = k)
    y_p_order <- lagm(as.matrix(y), p_order)
    colnames(y_p_order) <- gsub(x = colnames(y_p_order), 
                                pattern = "y", replacement = dep_var)
    xp1_q_order <- lagm(as.matrix(xp1), q_order1)
    xp2_q_order <- lagm(as.matrix(xp2), q_order2)
    
    xn1_q_order <- lagm(as.matrix(xn1), q_order1)
    xn2_q_order <- lagm(as.matrix(xn2), q_order2)
    
    {
      if (p_order > 1) {
        y_p_order <- y_p_order[-1, ]
      }
      else {
        y_p_order <- na.omit(y_p_order)
      }
    }
    y <- data.frame(y[-1])
    colnames(y) <- dep_var
    nardl_data <- data.frame(cbind(y, y_p_order, 
                                   xp1,xn1,  
                                   xp1_q_order, xn1_q_order,
                                   xp2, xn2, 
                                   xp2_q_order, xn2_q_order,
                                   trend))
    nardl_case <- names(nardl_data)
    lagx1 <- lagm(as.matrix(x1), 1)
    lagx2 <- lagm(as.matrix(x2), 1)
    
    lagx1 <- data.frame(lagx1[-1])
    lagx2 <- data.frame(lagx2[-1])
    
    colnames(lagx1) <- colnames(lagx1)
    colnames(lagx2) <- colnames(lagx2)
    
    dx1_qlag <- lagm(as.matrix(dx1), q_order1)
    dx2_qlag <- lagm(as.matrix(dx2), q_order2)
  }
  data <- rhs$data
  case_ <- names(data)
  
  {
    if (is.null(case) == FALSE){
      if (case == 1){
        nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(paste(nardl_case[-length(nardl_case)][-1], collapse = "+"), "-1")))
        ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(paste(case_[-length(case_)][-1], collapse = "+"), "-1")))
      }
      if (case == 3) {
        nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(nardl_case[-length(nardl_case)][-1], collapse = "+")))
        ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(case_[-length(case_)][-1], collapse = "+")))
      }
      if (case == 5){
        nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(nardl_case[-1], collapse = "+")))
        ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(case_[-1], collapse = "+")))
      }
      if (case == 2){
        nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(nardl_case[-length(nardl_case)][-1], collapse = "+")))
        ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(case_[-length(case_)][-1], collapse = "+")))
      }
      else if (case == 4){
        nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(nardl_case[-1], collapse = "+")))
        ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(case_[-1], collapse = "+")))
      }
    }
    if (is.null(case) == TRUE) {
      nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), paste(nardl_case[-1], collapse = "+")))
      ecm_fmla <- as.formula(paste(paste0("D.", dep_var, " ~ "), paste(case_[-1], collapse = "+")))
    }
  }
  
  lrname1 <- c(paste0(decomp1, "_pos_1"), paste0(decomp1, "_neg_1"))
  lrname2 <- c(paste0(decomp2, "_pos_1"), paste0(decomp2, "_neg_1"))
  
  {
    if(!is.null(case) == TRUE){
      if(gets == TRUE){
        stop("Case should be set as NULL when gets = TRUE. When case is an integer (1,2,3,4 or 5), gets should be set as FALSE")
      }
    }
    if(is.null(case) == TRUE){
      if (gets == FALSE)
        stop('Case should be an integer (1,2,3,4 or 5) when gets is set as FALSE')
    }
  }
  
  {
    if(gets == TRUE){
      nardl_fit_case_ <- lm(nardl_fmla, data = nardl_data, na.action = na.exclude)
      nardl_gets_fit <- gets.lm(nardl_fit_case_, wald.pval = gets_pval, include.1cut = T, print.searchinfo = F)
      
      fit_case_ <- lm(ecm_fmla, data = data, na.action = na.exclude)
      ecm_gets_fit <- gets.lm(fit_case_, wald.pval = gets_pval, 
                              keep = c(2:(1 + 1 + length(lrname1)+ length(lrname2))), include.1cut = T, print.searchinfo = F)
    }
    else{
      nardl_gets_fit <- lm(nardl_fmla, data = nardl_data, na.action = na.exclude)
      
      ecm_gets_fit <- lm(ecm_fmla, data = data, na.action = na.exclude)
    }
  }
  
  ecm_gets_summ <- summary(ecm_gets_fit)
  
  if (gets == TRUE){
    case_6 <- ifelse(c("(Intercept)", "trend") %in% rownames(ecm_gets_summ$coefficients), 1, 0)
    case = c()
    {
      if (sum(case_6) == 0) {
        case = 1
      }
      if (sum(case_6) == 2) {
        case = c(4, 5)
      }
      if (sum(case_6) == 1) {
        if (c("(Intercept)") %in% rownames(ecm_gets_summ$coefficients)) 
          case <- c(2, 3)
        else if (c("trend") %in% rownames(ecm_gets_summ$coefficients)) 
          stop("\n    You final model has a trend and no intercept. You should have a model with no trend and intercept; with intercept and no trend or both intercept and trend. \n    Consider adjusting the values for either the p_order, q_order or both in other to examine the various case - 1 to 5. Alternatively, you may consider adopting gets_nardl_uecm()")
      }
    }
  }
  
  k <- rhs$k
  summary_fit <- summary(ecm_gets_fit)
  coeff <- summary_fit$coefficients
  nlvars <- length(coeff[, 1])
  b_pos_neg <- c("_pos", "_neg")
  bp <- str_subset(rownames(coeff), b_pos_neg[1])[1:2]
  bn <- str_subset(rownames(coeff), b_pos_neg[2])[1:2]
  fstat_name = NULL
  {
    if (length(case) < 2) {
      if (case == 1 | case == 3 | case == 5) {
        fstat_name <- c(paste0(dep_var, "_1"), bn,  bp)
      }
      else if (case == 2) {
        fstat_name <- c("(Intercept)", paste0(dep_var, "_1"), bn, bp)
      }
      else {
        fstat_name <- c("trend", paste0(dep_var, "_1"), bn, bp)
      }
    }
    if (length(case) == 2) {
      for (i in 1:2) {
        if (case[[i]] == 1 | case[[i]] == 3 | case[[i]] == 5) {
          fstat_name[[i]] <- c(paste0(dep_var, "_1"), bn, bp)
        }
        else if (case[[i]] == 2) {
          fstat_name[[i]] <- c("(Intercept)", paste0(dep_var, "_1"), bn, bp) 
        }
        else {
          fstat_name[[i]] <- c("trend", paste0(dep_var, "_1"), bn, bp)
        }
      }
    }
  }
  
  {
    if (gets == TRUE){
      names(fstat_name) <- paste("Case", case)
    }
    else {fstat_name = fstat_name}
  }
  
  {
    if ("(Intercept)" %in% rownames(coeff) == TRUE) {
      lvars <- coeff[3:nlvars, 1]
      coof <- -lvars/coeff[[2]]
    }
    else {
      lvars <- coeff[2:nlvars, 1]
      coof <- -lvars/coeff[[1]]
    }
  }
  seldata <- data.matrix(coeff)
  cof <- matrix(coof, length(lvars), 1)
  
  {
    if ("(Intercept)" %in% rownames(coeff) == TRUE) {
      fb1 <- lvars/coeff[[2]]^2
      fb2 <- (-1/coeff[[2]]) * diag(nrow(as.matrix(fb1)))
    }
    else {
      fb1 <- lvars/coeff[[1]]^2
      fb2 <- (-1/coeff[[1]]) * diag(nrow(as.matrix(fb1)))
    }
  }
  fb <- cbind(as.matrix(fb1), fb2)
  vc <- vcov(ecm_gets_fit)
  
  {
    if ("(Intercept)" %in% rownames(coeff) == TRUE) {
      vcc <- vc[2:nrow(vc), 2:ncol(vc)]
    }
    else {
      vcc <- vc[1:nrow(vc), 1:ncol(vc)]
    }
  }
  lrse <- sqrt(diag(fb %*% vcc %*% t(fb)))
  lrt <- coof/lrse
  X <- sum(ifelse(names(ecm_gets_fit$model) %in% names(coef(ecm_gets_fit)), 1, 0))
  Y <- dim(ecm_gets_fit$model)[1]
  rdf <- Y - X
  lrpv <- 2 * pt(abs(coof/lrse), df = rdf, lower.tail = FALSE)
  (lres <- cbind(coof, lrse, lrt, lrpv))
  colnames(lres) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  (lres <- lres[c(bp, bn), ])
  (tstat <- coef(summary(ecm_gets_fit))[paste0(dep_var, "_1"), 3])
  wld_test <- c()
  fstat <- c()
  {
    if (length(fstat_name) > 2) {
      wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = fstat_name, 
                                   verbose = F, rhs = rep(0, length(fstat_name)), 
                                   test = "F")
      fstat <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2], 
                     df = wld_test$Df[2])
      rownames(fstat) <- c("wald_test")
      fstat
    }
    if (length(fstat_name) == 2) {
      for (i in 1:2) {
        wld_test[[i]] <- linearHypothesis(ecm_gets_fit, 
                                          hypothesis.matrix = fstat_name[[i]], verbose = F, 
                                          rhs = rep(0, length(fstat_name[[i]])), test = "F")
        fstat[[i]] <- cbind(Fstat = wld_test[[i]]$F[2], 
                            Pval = wld_test[[i]]$`Pr(>F)`[2])
        rownames(fstat[[i]]) <- c("wald_test")
      }
      names(fstat) <- paste("Case", case)
    }
  }
  ecm_diag_gets <- round(diag <- diagnostics(ecm_gets_summ, 
                                             ar.LjungB = c(floor(max((q_order1+q_order2)/2)/2), 0.025), arch.LjungB = NULL, 
                                             normality.JarqueB = NULL), 3)
  jbtest <- c()
  jbtest$statistic <- round(jarque.bera.test(ecm_gets_fit$residuals)$statistic, 4)
  jbtest$p.value <- round(jarque.bera.test(ecm_gets_fit$residuals)$p.value, 4)
  jbtest <- cbind(jbtest)
  jbtest <- cbind(jbtest[[1]], jbtest[[2]])
  colnames(jbtest) <- c("statistics", "p.value")
  rownames(jbtest) <- dep_var
  lm_test <- c()
  lm_test$statistic <- round(bgtest(ecm_gets_fit, type = "Chisq", 
                                    order = order_l)$statistic, 4)
  lm_test$p.value <- round(bgtest(ecm_gets_fit, type = "Chisq", 
                                  order = order_l)$p.value, 4)
  lm_test <- cbind(lm_test)
  lm_test <- cbind(lm_test[[1]], lm_test[[2]])
  colnames(lm_test) <- c("statistics", "p.value")
  rownames(lm_test) <- dep_var
  arch <- c()
  arch$statistic <- round(ArchTest(ecm_gets_fit$residuals, 
                                   order_l)$statistic, 4)
  arch$p.value <- round(ArchTest(ecm_gets_fit$residuals, order_l)$p.value, 4)
  arch <- cbind(arch)
  arch <- cbind(arch[[1]], arch[[2]])
  colnames(arch) <- c("statistics", "p.value")
  rownames(arch) <- dep_var
  reset_test <- c()
  reset_test$statistic <- round(resettest(ecm_gets_fit, power = 2, 
                                          type = "princomp")$statistic, 4)
  reset_test$p.value <- round(resettest(ecm_gets_fit, power = 2, 
                                        type = "princomp")$p.value, 4)
  reset_test <- cbind(reset_test)
  reset_test <- cbind(reset_test[[1]], reset_test[[2]])
  colnames(reset_test) <- c("statistics", "p.value")
  rownames(reset_test) <- dep_var
  diag <- rbind(lm_test, arch, jbtest, reset_test)
  rownames(diag) <- c("BG_SC_lm_test", "LM_ARCH_test", "normality_test", "RESET_test")
  coof_lres <- coof[c(bn, bp)]
  
  vcc[c(names(coof_lres[1:2]), names(coof_lres[3:4])), c(names(coof_lres[1:2]), names(coof_lres[3:4]))]
  
  vcc1 <- vcc[c(names(coof_lres[1:2]), names(coof_lres[3:4])), 
              c(names(coof_lres[1:2]), names(coof_lres[3:4]))]
  
  lbeta <- coeff[c(names(coof_lres[1:2]), names(coof_lres[3:4])), 1]
  ll2 <- c(paste(names(coof_lres[1:2]), "=", names(coof_lres[3:4])))
  wld_test1 <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = ll2[1], verbose = F, test = "F")
  lr_asym_test1 <- cbind(Fstat = wld_test1$F[2], Pval = wld_test1$`Pr(>F)`[2])
  
  wld_test2 <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = ll2[2], verbose = F, test = "F")
  lr_asym_test2 <- cbind(Fstat = wld_test2$F[2], Pval = wld_test2$`Pr(>F)`[2])
  rownames(lr_asym_test1) <- decomp1
  rownames(lr_asym_test2) <- decomp2
  
  decomp_d1 = paste0("D.", decomp1)
  decomp_d2 = paste0("D.", decomp2)
  
  sr_cof_nm <- str_subset(rownames(coeff), decomp_d1)
  coeff_sr <- coeff[sr_cof_nm, ]
  don <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_neg")))[1]
  dop <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_pos")))[1]
  
  {
    if (gets == TRUE){
      nnn <- str_subset(sr_cof_nm, pattern = "_neg")
      ppp <- str_subset(sr_cof_nm, pattern = "_pos")
    }
    else { 
      nnn = NULL
      ppp = NULL
    }
  }
  
  neg <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_neg"), " "))
  pos <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_pos"), " "))
  
  {
    if (don > 0 & dop > 0) {
      neg = paste(str_split(neg, boundary("word"))[[1]], "+")
      pos = paste(str_split(pos, boundary("word"))[[1]], "+")
      pos = c("=", pos)
      pos = c(pos[-length(pos)], str_split(pos[length(pos)], 
                                           boundary("word"))[[1]])
      neg = c(neg[-length(neg)], str_split(neg[length(neg)], 
                                           boundary("word"))[[1]])
      pos_neg <- str_flatten(c(neg, pos))
      wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_neg, 
                                   rhs=NULL, verbose = F, test = "F")
      sr_asym_test1 <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
      rownames(sr_asym_test1) <- decomp1
    }
    else {
      if (length(ppp) > 0) {
        message("Only one of the differenced - decomposed variable (decomp1_) appeared in the final model. The value in `Shortrun_asym` is the wald test on the sum of the coefficients of the available short-run decomposed variable does not have any significant effect on the best model")
        pos = paste(str_split(pos, boundary("word"))[[1]], 
                    "+")
        pos = c(pos[-length(pos)], str_split(pos[length(pos)], 
                                             boundary("word"))[[1]])
        pos = c(pos, "= 0")
        pos_h <- str_flatten(pos)
        wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_h, 
                                     rhs=NULL, verbose = F, test = "F")
        sr_asym_test1 <- cbind(Fstat = wld_test$F[2], 
                               Pval = wld_test$`Pr(>F)`[2])
        rownames(sr_asym_test1) <- decomp1
      }
      else {
        if (length(nnn) > 0) {
          neg = paste(str_split(neg, boundary("word"))[[1]], 
                      "+")
          neg = c(neg[-length(neg)], str_split(neg[length(neg)], 
                                               boundary("word"))[[1]])
          neg = c(neg, "= 0")
          neg_h <- str_flatten(neg)
          wld_test <- linearHypothesis(ecm_gets_fit, 
                                       hypothesis.matrix = neg_h, rhs=NULL, verbose = F, 
                                       test = "F")
          message("Only one of the differenced - decomposed variable (decomp1_) appeared in the final model. The value in `Shortrun_asym` is the wald test on the sum of the coefficients of the available short-run decomposed variable does not have any significant effect on the best model")
          sr_asym_test1 <- cbind(Fstat = wld_test$F[2], 
                                 Pval = wld_test$`Pr(>F)`[2])
          rownames(sr_asym_test1) <- decomp1
        }
        else {
          message(("No short-run coefficients for the first decomposed variable "))
          sr_asym_test1 = c(("No short-run coefficients for the second decomposed variable"))
        }
      }
    }
  }
  
  sr_cof_nm <- str_subset(rownames(coeff), decomp_d2)
  coeff_sr <- coeff[sr_cof_nm, ]
  don <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_neg")))[1]
  dop <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_pos")))[1]
  
  {
    if (gets == TRUE){
      nnn <- str_subset(sr_cof_nm, pattern = "_neg")
      ppp <- str_subset(sr_cof_nm, pattern = "_pos")
    }
    else { 
      nnn = NULL
      ppp = NULL
    }
  }
  
  neg <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_neg"), " "))
  pos <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_pos"), " "))
  
  {
    if (don > 0 & dop > 0) {
      (neg = paste(str_split(neg, boundary("word"))[[1]], "+"))
      (pos = paste(str_split(pos, boundary("word"))[[1]], "+"))
      (pos = c("=", pos))
      pos = c(pos[-length(pos)], str_split(pos[length(pos)],boundary("word"))[[1]])
      neg = c(neg[-length(neg)], str_split(neg[length(neg)], boundary("word"))[[1]])
      pos_neg <- str_flatten(c(neg, pos))
      wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_neg, verbose = F, test = "F")
      sr_asym_test2 <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
      rownames(sr_asym_test2) <- decomp2
    }
    else {
      if (length(ppp) > 0) {
        message("Only one of the differenced - decomposed variable (D.decomp2_pos) appears in the final model. The value in `Shortrun_asym` is the wald test on the sum of the coefficients of the available short-run decomposed variable does not have any significant effect in the best model")
        pos = paste(str_split(pos, boundary("word"))[[1]], "+")
        pos = c(pos[-length(pos)], str_split(pos[length(pos)], boundary("word"))[[1]])
        pos = c(pos, "= 0")
        pos_h <- str_flatten(pos)
        wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_h, 
                                     verbose = F, test = "F")
        sr_asym_test2 <- cbind(Fstat = wld_test$F[2], 
                               Pval = wld_test$`Pr(>F)`[2])
        rownames(sr_asym_test2) <- decomp2
      }
      else {
        if (length(nnn) > 0) {
          neg = paste(str_split(neg, boundary("word"))[[1]], 
                      "+")
          neg = c(neg[-length(neg)], str_split(neg[length(neg)], 
                                               boundary("word"))[[1]])
          neg = c(neg, "= 0")
          neg_h <- str_flatten(neg)
          wld_test <- linearHypothesis(ecm_gets_fit, 
                                       hypothesis.matrix = neg_h, verbose = F, 
                                       test = "F")
          message("Only one of the differenced - decomposed variable (D.decomp2_neg) appears in the final model. The value in `Shortrun_asym` is the wald test on the sum of the coefficients of the available short-run decomposed variable does not have any significant effect on the best model")
          sr_asym_test2 <- cbind(Fstat = wld_test$F[2], 
                                 Pval = wld_test$`Pr(>F)`[2])
          rownames(sr_asym_test2) <- decomp2
        }
        else {
          message("No short-run coefficients for the second decomposed variable")
          sr_asym_test2 = c("No short-run coefficients for the second decomposed variable")
        }
      }
    }
  }
  
  nobs <- nobs(ecm_gets_fit)
  e <- ecm_gets_fit$residuals
  stab_plot <- function(graph_save) {
    if (graph_save == TRUE) {
      oldpar <- par(no.readonly = TRUE)
      e <- ecm_gets_fit$residuals
      n <- nobs
      par(mfrow = c(1, 2))
      cusum(e = e, k = k, n = n)
      cumsq(e = e, k = k, n = n)
      on.exit(par(oldpar))
    }
  }
  stab_plot(graph_save)
  
  bounds_F <- c()
  {
    if (length(case) < 2) {
      {
        if (case == 1 | case == 3 | case == 5) {
          bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                              fstat = fstat[1, 1], tstat = tstat, obs = nobs, 
                                              k = k)
        }
        else if (case == 2) {
          bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                              fstat = fstat[1, 1], tstat = NULL, obs = nobs, 
                                              k = k)
        }
        else {
          bounds_F <- dynamac_pkg_bounds_test(case = case, 
                                              fstat = fstat[1, 1], tstat = NULL, obs = nobs, 
                                              k = k)
        }
      }
    }
    if (length(case) == 2) {
      for (i in 1:2) {
        if (case[[i]] == 1 | case[[i]] == 3 | case[[i]] == 
            5) {
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case = case[[i]], 
                                                   fstat = fstat[[i]][1, 1], tstat = tstat, 
                                                   obs = nobs, k = k)
        }
        else if (case[[i]] == 2) {
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case = case[[i]], 
                                                   fstat = fstat[[i]][1, 1], tstat = NULL, 
                                                   obs = nobs, k = k)
        }
        else {
          bounds_F[[i]] <- dynamac_pkg_bounds_test(case = case[[i]], 
                                                   fstat = fstat[[i]][1, 1], tstat = NULL, 
                                                   obs = nobs, k = k)
        }
      }
      names(bounds_F) <- paste("Case", case)
    }
  }
  lr_asym_test <- rbind(lr_asym_test1, lr_asym_test2)
  sr_asym_test <- list(sr_asym_test1,sr_asym_test2)
  names(sr_asym_test) <- c(decomp1,decomp2)
  gets_NARDL <- list(NARDL_fit = nardl_gets_fit, 
                     ECM_fit = ecm_gets_fit, 
                     Summary_uecm_fit = ecm_gets_summ, 
                     ecm_diagnostics_test = diag, 
                     longrun_asym = lr_asym_test,
                     Shortrun_asym = sr_asym_test, 
                     cointegration = bounds_F, 
                     Longrun_relation = lres)
  class(gets_NARDL) <- 'NARDL-gets'
  return(gets_NARDL)
}
