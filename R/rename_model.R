output_ren <- function(x,listn=listn,D.patern,D.repl,l.patern,l.repl)
  {
  lr_names <-lapply(1:length(listn), function(i) rownames(x[[i]]$Longrun_relation))
  lr_newnames_sr <- c()
  lr_newnames_sr <- lapply(1:length(lr_names), function(i)
    gsub(x = lr_names,  pattern = D.patern,  replacement = D.repl))
  lr_newnames_sr <- lapply(1:length(lr_names), function(i)
    gsub(x = lr_names[[i]],  pattern = l.patern,  replacement = l.repl))
  names(lr_newnames_sr) <- listn
  
  for (i in 1:length(lr_names)) {
    rownames(x[[i]]$Longrun_relation) <- c(lr_newnames_sr[[i]])
  }
  uecm_names <-lapply(1:length(listn), function(i) 
    rownames(x[[i]]$UECM$coefficients))
  
  uecm_new_names <- lapply(1:length(uecm_names), function(i)
    gsub(x = uecm_names,  pattern = D.patern,  replacement = D.repl))
  uecm_new_names <- lapply(1:length(uecm_names), function(i)
    gsub(x = uecm_names[[i]],  pattern = l.patern,  replacement = l.repl))
  names(uecm_new_names) <- listn
  for (i in 1:length(uecm_names)) {
    rownames(x[[i]]$UECM$coefficients) <- c(uecm_new_names[[i]])
  }
  return(x)
}