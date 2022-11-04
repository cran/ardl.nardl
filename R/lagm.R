lagm <- function(x,k){
  dc <- c()
  {
    if(k >= 1){
      if(dim(x)[2] > 1)
        for (j in 1:dim(x)[2]) {
          dc[[j]] <- lapply(1:k, function(i) lag(x[,j],i))
        }
      if(dim(x)[2] == 1){
        dc <- lapply(1:k, function(i) lag(x,i))
      }
    }
    if(k == 0){
      dc <- x
    }
    if(dim(x)[2] > 1){
      (dc <-lapply(1:dim(x)[2], function(i) list.cbind(dc[[i]])))
      (dc <- list.cbind(dc))
      colnames(dc) <- unlist(lapply(1:dim(x)[2], function(i) paste0(colnames(x)[i],'_',1:k)))
    }
    else{
      (dc <- list.cbind(dc))
      colnames(dc) <- paste0(colnames(x),'_',1:k)
    }
  }
  return(dc)
}
