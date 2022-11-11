deduceHardCallFromGenotypeProbabilities <- function(aa, ab, bb, max.uncertainty  = NULL){
  if(!is.null(max.uncertainty) && 1 - max(aa, ab, bb) > max.uncertainty){
    return(NA)
  }
  
  CALLS <- c(0, 1, 2)
  return(CALLS[which.max(c(aa, ab, bb))])
}

collapseGenotypeTriplets <- function(x, max.uncertainty  = NULL){
  stopifnot(length(x) %% 3 == 0)
  rs <- sapply(
    seq(1, length(x), 3),
    function(first.idx){
      probs <- as.list(x[first.idx:(first.idx + 2)])
      names(probs) <- c("aa", "ab", "bb")
      probs[["max.uncertainty"]] <- max.uncertainty
      return(do.call(deduceHardCallFromGenotypeProbabilities, probs))
    }
  )
  
  return(rs)
}

collapseGenotypeTripletMatrix <- function(x, max.uncertainty = NULL){
  rs <- apply(
    x,
    1,
    collapseGenotypeTriplets,
    max.uncertainty
  )
  
  return(t(rs))
}

oxfordGenToHardCalls <- function(x, max.uncertainty = NULL, chr = NULL){
  last.meta.col.idx <- 5
  if((ncol(x) - last.meta.col.idx) %% 3 != 0){
    last.meta.col.idx <- 6
  }
  
  stopifnot((ncol(x) - last.meta.col.idx) %% 3 == 0)
  df.meta <- x[, 1:last.meta.col.idx]
  df.probs <- x[, (last.meta.col.idx + 1):ncol(x)]
  df.probs <- collapseGenotypeTripletMatrix(df.probs, max.uncertainty)
  stopifnot(nrow(df.meta) == nrow(df.probs))
  rs <- cbind(df.meta, df.probs)
  if(!is.null(chr)){
    rs[, 1] <- chr
  }
  
  return(rs)
}




#####
# QUICK TEST
#####
(function(){
  library(dplyr)
  test.data <- "
  1	rs1231	100	A	G	0.2	0.8	0	0.3	0.1	0.6	0.4	0.4	0.2
  4	rs1232	101	T	C	0.9	0.05	0.05	0.4	0.55	0.05	0.3	0.1	0.6
  2	rs1233	102	G	A	0	0.4	0.6	0.8	0.1	0.1	0.1	0.8	0.1
  1	rs1234	103	C	T	0.95	0.04	0.01	0.1	0.9	0	0	0.2	0.8
  "
  test.data <- paste(sub("^\\s*", "", readLines(textConnection(test.data))), collapse = "\n")
  
  expected <- "
  1	rs1231	100	A	G	1	2	0
  4	rs1232	101	T	C	0	1	2
  2	rs1233	102	G	A	2	0	1
  1	rs1234	103	C	T	0	1	2
  "
  df.expected = data.frame(
    V1 = c(1,4,2,1),
    V2 = c("rs1231", "rs1232", "rs1233", "rs1234"),
    V3 = c(100, 101, 102, 103),
    V4 = c("A", "T", "G", "C"),
    V5 = c("G", "C", "A", "T"),
    `1` = c(1, 0, 2, 0),
    `2` = c(2, 1, 0, 1),
    `3` = c(0, 2, 1, 2),
    check.names = F
  )
  
  
  df.test.data <- read.table(textConnection(test.data), sep = "\t")
  df.actual <- oxfordGenToHardCalls(df.test.data, NULL)
  stopifnot(all.equal(df.actual, df.expected))
  print("SUCCESS")
})()
