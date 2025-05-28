# evaluate the target quantities at a range of lower thresholds

Q1diag <- function(sims) {
  x <- seq(1,1.7,by=0.1)
#  tvc <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims1,FUN = Q1,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
#  temp <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims2,FUN = Q1,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
  tvc <- do.call(rbind,lapply(sims1,FUN = Q1,threshold=x[i]))

    # calculate empirical values
  tmp <- rbind(tvc,temp)
  colnames(tmp) <- as.character(x)
  return(tmp) 
}

Q2diag <- function(sims1,sims2,Nrun=50) {
  x <- seq(2.1,5.7,by=0.6)
  tvc <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims1,FUN = Q2,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
  temp <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims2,FUN = Q2,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
  # calculate empirical values
  tmp <- rbind(tvc,temp)
  colnames(tmp) <- as.character(x)
  return(tmp) 
}

Q3diag <- function(sims1,sims2,Nrun=50) {
  x <- seq(2,5,by=0.5)
  tvc <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims1,FUN = Q3,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
  temp <- sapply(1:length(x),FUN= function(i) {apply(do.call(rbind,lapply(sims2,FUN = Q3,threshold=x[i])),MARGIN=2,FUN=sum)/Nrun})
  # calculate empirical values
  tmp <- rbind(tvc,temp)
  colnames(tmp) <- as.character(x)
  return(tmp)  
}

