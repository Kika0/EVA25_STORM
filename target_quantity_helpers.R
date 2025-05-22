Q1prel <- function(Yrun) {
  sum(apply(Yrun,MARGIN=1,FUN=sum)>85)
}

Q2prel <- function(Yrun) {
  sum(apply(X=apply(X=Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))),MARGIN = c(1,2),function(x){x>4.3}),MARGIN=1,sum)>=3)
}

Q3prel <- function(Yrun) {
  str_count(paste(c(0,0,apply(X=apply(X=Yrun %>% dplyr::select(all_of(c(1,6,11,16,21))),MARGIN = c(1,2),function(x){x>2.5}),MARGIN=1,sum)>=3,0,0),collapse=""),pattern="110")
}

Q1 <- function(Yrun,threshold=1.7) {
  sum(apply(apply(X=Yrun,MARGIN = c(1,2),FUN = function(x){x>threshold}),MARGIN=1,FUN=sum)==25)
}

Q2 <- function(Yrun,threshold=5.7) {
  sum(apply(X=apply(X=Yrun,MARGIN = c(1,2),function(x){x>threshold}),MARGIN=1,sum)>=6)
}

Q3 <- function(Yrun,threshold=5) {
  str_count(paste(c(0,0,apply(X=apply(X=Yrun,MARGIN = c(1,2),function(x){x>threshold}),MARGIN=1,sum)>=3,0,0),collapse=""),pattern="110")
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
