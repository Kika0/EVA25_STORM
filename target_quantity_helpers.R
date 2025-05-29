Q1prel <- function(Yrun, threshold = 85) {
  sum(rowSums(Yrun) > threshold)
}

# Q2prel <- function(Yrun) {
#   sum(apply(X=apply(X=Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))),MARGIN = c(1,2),function(x){x>4.3}),MARGIN=1,sum)>=3)
# }

Q2prel <- function(Yrun, threshold = 4.3) {
  sum(rowSums(Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))) > threshold) >= 3)
  
}

# Q3prel <- function(Yrun) {
#   str_count(paste(c(0,0,apply(X=apply(X=Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))),MARGIN = c(1,2),function(x){x>2.5}),MARGIN=1,sum)>=3,0,0),collapse=""),pattern="110")
# }

Q3prel <- function(Yrun, threshold = 2.5){
  
  ## Obtain the excesses for the chosen site
  excesses <- rowSums(Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))) > threshold) >= 3
  
  ## Turn the above into a string of 0's and 1's
  excesses_string <- paste(as.numeric(excesses), collapse = "")
  
  ## Count the number of times a run of two appears
  ## Method will only count run of more than two 1s once
  out <- str_count(excesses_string, "1{2,}")
  
  ## return the output
  return(out)
}

Q1 <- function(Yrun,threshold=1.7) {
  sum(rowSums(Yrun > threshold) == 25)
}

Q2 <- function(Yrun,threshold=5.7) {
  sum(rowSums(Yrun > threshold) >= 6)
}

Q3 <- function(Yrun, threshold = 5){
  ## Number of times at least three sites exceed the threshold
  excesses <- rowSums(Yrun > threshold) >= 3
  
  ## Turn the above into a string of 0's and 1's
  excesses_string <- paste(as.numeric(excesses), collapse = "")
  
  ## Count the number of times a run of two appears
  ## Method will only count run of more than two 1s once
  out <- str_count(excesses_string, "1{2,}")
  
  ## return the output
  return(out)
}

# evaluate all quantities from one run
Qeval <- function(sims) {
  return(data.frame("Q1prel"=Q1prel(sims),
                    "Q2prel"=Q2prel(sims),
                    "Q3prel"=Q3prel(sims),
                    "Q1"=Q1(sims),
                    "Q2"=Q2(sims),
                    "Q3"=Q3(sims)))
}
