gets_nardl_uecm <- function(x, 
                            decomp, 
                            dep_var, 
                            control = NULL, 
                            c_q_order = c(2), 
                            p_order = c(3), 
                            q_order = c(4),
                            gets_pval = 0.1, 
                            order_l = 4,
                            F_HC = FALSE,
                            graph_save = FALSE, 
                            case = 3){
{
  if (!is.null(control) == T) {
    x <- na.omit(x)
    df <- c()
    df$y <- x[, dep_var]
    df$x1 <- x[, decomp]
    df$x2 <- x[, control]
    x <- model.matrix(~x1 - 1, df)
    y <- model.matrix(~y - 1, df)
    colnames(x) <- decomp
    dy <- diff(y)
    colnames(dy) <- paste("D", dep_var, sep = ".")
    dx <- diff(x)
    colnames(dx) <- paste("D", colnames(x), sep = ".")
    h <- model.matrix(~x2 - 1, df)
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
    colnames(lagy) <- paste0(dep_var, "_1")
    lagh <- lagm(as.matrix(h), 1)
    colnames(lagh) <- paste0(control, "_1")
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
    rhs <- cbind(dy, lagy, lxp, lxn, lagh, ldy, ldxp, ldxn, 
                 ldh, trend)
    colnames(rhs)
    rhs <- list(data = data.frame(rhs), k = k)
    y_p_order <- lagm(as.matrix(y), p_order)
    colnames(y_p_order) <- gsub(x = colnames(y_p_order), 
                                pattern = "y", replacement = dep_var)
    xp_q_order <- lagm(as.matrix(xp), q_order)
    xn_q_order <- lagm(as.matrix(xn), q_order)
    h_c_q_order <- lagm(as.matrix(h), c_q_order)
    colnames(h_c_q_order) <- gsub(x = colnames(h_c_q_order), 
                                  pattern = "x2", replacement = control)
    {
      if (p_order > 1) {
        y_p_order <- y_p_order[-1, ]
      }
      else {
        y_p_order <- na.omit(y_p_order)
      }
    }
    {
      if (c_q_order > 1) {
        h_c_q_order <- h_c_q_order[-1, ]
      }
      else {
        h_c_q_order <- na.omit(h_c_q_order)
      }
    }
    y <- data.frame(y[-1])
    colnames(y) <- dep_var
    h <- data.frame(h[-1])
    colnames(h) <- control
    nardl_data <- data.frame(cbind(y, y_p_order, xp, xp_q_order, 
                                   xn, xn_q_order, h, h_c_q_order, trend))
    nardl_case <- names(nardl_data)
    lagx <- lagm(as.matrix(x), 1)
    lagx1 <- data.frame(lagx[-1])
    colnames(lagx1) <- colnames(lagx)
    rhs_lng_sym <- cbind(dy, lagy, lagx1, lagh, ldy, ldxp, 
                         ldxn, ldh, trend)
    dx_qlag <- lagm(as.matrix(dx), q_order)
    rhs_shr_sym <- cbind.data.frame(dy, lagy, lxp, lxn, 
                                    lagh, ldy, dx_qlag, ldh, trend)
  }
  else {
    x <- na.omit(x)
    df <- c()
    df$y <- x[, dep_var]
    df$x1 <- x[, decomp]
    x <- model.matrix(~x1 - 1, df)
    y <- model.matrix(~y - 1, df)
    colnames(x) <- decomp
    colnames(y) <- dep_var
    lhs <- "y"
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
    rhs <- cbind(dy, lagy, lxp, lxn, ldy, ldxp, ldxn, trend)
    rhs <- list(data = data.frame(rhs), k = k)
    y_p_order <- lagm(as.matrix(y), p_order)
    colnames(y_p_order) <- gsub(x = colnames(y_p_order), 
                                pattern = "y", replacement = dep_var)
    xp_q_order <- lagm(as.matrix(xp), q_order)
    xn_q_order <- lagm(as.matrix(xn), q_order)
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
    nardl_data <- data.frame(cbind(y, y_p_order, xp, xp_q_order, 
                                   xn, xn_q_order, trend))
    nardl_case <- names(nardl_data)
    lagx <- lagm(as.matrix(x), 1)
    lagx1 <- data.frame(lagx[-1])
    colnames(lagx1) <- colnames(lagx)
    rhs_lng_sym <- cbind(dy, lagy, lagx1, ldy, ldxp, ldxn, 
                         trend)
    dx_qlag <- lagm(as.matrix(dx), q_order)
    rhs_shr_sym <- cbind.data.frame(dy, lagy, lxp, lxn, 
                                    ldy, dx_qlag, trend)
  }
}
data <- rhs$data
case_ <- names(data)
{
  if (case == 1) 
    nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), 
                                   paste(paste(nardl_case[-length(nardl_case)][-1], 
                                               collapse = "+"), "-1")))
  if (case == 3) 
    nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), 
                                   paste(nardl_case[-length(nardl_case)][-1], collapse = "+")))
  if (case == 5) 
    nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), 
                                   paste(nardl_case[-1], collapse = "+")))
  if (case == 2) 
    nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), 
                                   paste(nardl_case[-length(nardl_case)][-1], collapse = "+")))
  if (case == 4) 
    nardl_fmla <- as.formula(paste(paste0(dep_var, " ~ "), 
                                   paste(nardl_case[-1], collapse = "+")))
}
{
  if (case == 1) 
    ecm_fmla <- as.formula(paste(paste0("D.", dep_var, 
                                        " ~ "), paste(paste(case_[-length(case_)][-1], 
                                                            collapse = "+"), "-1")))
  if (case == 3) 
    ecm_fmla <- as.formula(paste(paste0("D.", dep_var, 
                                        " ~ "), paste(case_[-length(case_)][-1], collapse = "+")))
  if (case == 5) 
    ecm_fmla <- as.formula(paste(paste0("D.", dep_var, 
                                        " ~ "), paste(case_[-1], collapse = "+")))
  if (case == 2) 
    ecm_fmla <- as.formula(paste(paste0("D.", dep_var, 
                                        " ~ "), paste(case_[-length(case_)][-1], collapse = "+")))
  if (case == 4) 
    ecm_fmla <- as.formula(paste(paste0("D.", dep_var, 
                                        " ~ "), paste(case_[-1], collapse = "+")))
}

nardl_fit_case_ <- lm(nardl_fmla, data = nardl_data, na.action = na.exclude)

nardl_gets_fit <- gets.lm(nardl_fit_case_,
                          wald.pval=gets_pval,
                          include.1cut = T,print.searchinfo = F)


fit_case_ <- lm(ecm_fmla, data = data, na.action = na.exclude)

{
  if(is.null(control) == FALSE){
    (lrname <- c(paste0(decomp,'_pos_1'), paste0(control,'_1'), paste0(decomp,'_neg_1')))
  }
  else{
    (lrname <- c(paste0(decomp,'_pos_1'), paste0(decomp,'_neg_1')))
  }
}

{
  if("(Intercept)" %in% names(fit_case_$coefficients) == TRUE){
    keep_val <- c(1:(1+1+length(lrname)))
  }
  else{
    keep_val <- c(1:(1+length(lrname)))
  }
}

{
  if("trend" %in% names(fit_case_$coefficients) == TRUE){
    keep_val <- c(keep_val,length(fit_case_$coefficients))
  }
  else{
    keep_val <- keep_val
  }
}

ecm_gets_fit <- gets.lm(fit_case_,
                        wald.pval=gets_pval,
                        keep = keep_val,
                        include.1cut = T,print.searchinfo = F)
k <- rhs$k
summary_fit <- summary(ecm_gets_fit)
coeff <- summary_fit$coefficients
nlvars <- length(coeff[, 1])
b_pos_neg <- c("_pos", "_neg")
bp <- str_subset(rownames(coeff), b_pos_neg[1])[1]
bn <- str_subset(rownames(coeff), b_pos_neg[2])[1]
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
X <- sum(ifelse(names(ecm_gets_fit$model) %in% names(coef(ecm_gets_fit)), 
                1, 0))
Y <- dim(ecm_gets_fit$model)[1]
rdf <- Y - X
lrpv <- 2 * pt(abs(coof/lrse), df = rdf, lower.tail = FALSE)
lres <- cbind(coof, lrse, lrt, lrpv)
colnames(lres) <- c("Estimate", "Std. Error", "t value", 
                    "Pr(>|t|)")
{
  if (is.null(control) == FALSE) {
    control_lag_v <- paste0(control, "_1")
    lres <- lres[c(bp, bn, control_lag_v), ]
  }
  else {
    lres <- lres[c(bp, bn), ]
  }
}
tstat <- coef(summary(ecm_gets_fit))[paste0(dep_var, "_1"), 
                                     3]
{
  if (!is.null(control) == T) {
    if (case == 1 | case == 3 | case == 5) {
      fstat_name <- c(paste0(dep_var, "_1"), bn, bp, 
                      paste0(control, "_1"))
    }
    else if (case == 2) {
      fstat_name <- c("(Intercept)", paste0(dep_var, 
                                            "_1"), bn, bp, paste0(control, "_1"))
    }
    else {
      fstat_name <- c("trend", paste0(dep_var, "_1"), 
                      bn, bp, paste0(control, "_1"))
    }
  }
  else {
    if (case == 1 | case == 3 | case == 5) {
      fstat_name <- c(paste0(dep_var, "_1"), bn, bp)
    }
    else if (case == 2) {
      fstat_name <- c("(Intercept)", paste0(dep_var, 
                                            "_1"), bn, bp)
    }
    else {
      fstat_name <- c("trend", paste0(dep_var, "_1"), 
                      bn, bp)
    }
  }
}

wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = fstat_name, verbose = F,
                              rhs = rep(0,length(fstat_name)), test = "F")
wld_test_ <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = fstat_name, verbose = F,
                               rhs = rep(0,length(fstat_name)), test = "F",vcov = vcovHC)
lmwld <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2], df = wld_test$Df[2])
lmwld_ <- cbind(Fstat = wld_test_$F[2], Pval = wld_test_$`Pr(>F)`[2], df = wld_test_$Df[2])

fstat <- rbind(lmwld,lmwld_)
rownames(fstat) <- c('wald_test','robust_wald_test')
ecm_gets_summ <- summary(ecm_gets_fit)
ecm_diag_gets <- round(diag <- diagnostics(ecm_gets_summ,
                                            ar.LjungB=c(floor(max(q_order)/2), 0.025),
                                            arch.LjungB = NULL, 
                                            normality.JarqueB = NULL),4)

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

coof_lres <- coof[c(bn, bp)]
vcc1 <- vcc[c(names(coof_lres[1]), names(coof_lres[2])), 
            c(names(coof_lres[1]), names(coof_lres[2]))]
lbeta <- coeff[c(names(coof_lres[1]), names(coof_lres[2])), 1]

ll2 <- c(paste(names(coof_lres[1]), "=", names(coof_lres[2])))

wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = ll2, verbose = F,test = "F")
lr_asym_test <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
rownames(lr_asym_test) <- decomp

decomp_d = paste0("D.", decomp)
sr_cof_nm <- str_subset(rownames(coeff), decomp_d)
coeff_sr <- coeff[sr_cof_nm, ]

don <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_neg")))[1]
dop <- dim(as.data.frame(str_subset(sr_cof_nm, pattern = "_pos")))[1]
nnn <- str_subset(sr_cof_nm, pattern = "_neg")
ppp <- str_subset(sr_cof_nm, pattern = "_pos")

neg <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_neg"), ' '))
pos <- str_flatten(paste(str_subset(sr_cof_nm, pattern = "_pos"), ' '))
{
  if (don > 0 & dop > 0) {
    neg = paste(str_split(neg, boundary("word"))[[1]], "+")
    pos = paste(str_split(pos, boundary("word"))[[1]], "+")
    pos = c("=", pos)
    pos = c(pos[-length(pos)], str_split(pos[length(pos)], boundary("word"))[[1]])
    neg = c(neg[-length(neg)], str_split(neg[length(neg)], boundary("word"))[[1]])
    pos_neg <- str_flatten(c(neg, pos))
    wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_neg, verbose = F, test = "F")
    sr_asym_test <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
    rownames(sr_asym_test) <- decomp
  }
  else {
    if (length(ppp) > 0) {
      pos = paste(str_split(pos, boundary("word"))[[1]], "+")
      pos = c(pos[-length(pos)], str_split(pos[length(pos)], boundary("word"))[[1]])
      pos = c(pos, "= 0")
      pos_h <- str_flatten(pos)
      wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = pos_h, verbose = F, test = "F")
      sr_asym_test <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
      rownames(sr_asym_test) <- decomp
    }
    else {
      if (length(nnn) > 0) {
        neg = paste(str_split(neg, boundary("word"))[[1]], "+")
        neg = c(neg[-length(neg)], str_split(neg[length(neg)], boundary("word"))[[1]])
        neg = c(neg, "= 0")
        neg_h <- str_flatten(neg)
        wld_test <- linearHypothesis(ecm_gets_fit, hypothesis.matrix = neg_h, verbose = F, test = "F")
        sr_asym_test <- cbind(Fstat = wld_test$F[2], Pval = wld_test$`Pr(>F)`[2])
        rownames(sr_asym_test) <- decomp
      }
      else {
        sr_asym_test = c('This model is similar to Short-run symmetric restriction (SRSR). Thus, no need for short-run asymmetric test. See nardl_uecm_sym() for more details.')
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
gets_NARDL <- list('Parsimonious_NARDL_fit' = nardl_gets_fit, 
                   'Parsimonious_ECM_fit' = ecm_gets_fit, 
                   'Summary_uecm_fit' = ecm_gets_summ,
                   'ecm_diagnostics_test' = diag,
                   'longrun_asym' = lr_asym_test, 
                   'Shortrun_asym' = sr_asym_test,
                   'cointegration'= bounds_F,
                   'Longrun_relation' = lres)
return(gets_NARDL)
}
