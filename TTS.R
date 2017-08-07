TTS <- function(x, ...) {
  UseMethod("TTS", x)
}

TTS.numeric <- function(stregnth_dist_properties = c(1, .5),
                        list_length = 10,
                        threshold = 0,
                        duration = 30,
                        constant_rate = TRUE) {
  library(whoppeR)
  print("numeric")
}

TTS.data.frame <- function(.data, ...) {
  
  is_column_numeric <- vapply(.data, is.numeric, FUN.VALUE = logical(1))
  if (!all(is_column_numeric)) {
    stop("All variables must be numeric")
  }

  x <- TTS.matrix(as.data.frame(.data), ...)
  return(x)
}

TTS.matrix <- function(.data,
                       duration = 30,
                       constant_rate = TRUE) {
  
  sim_RT  <- matrix(NA,
                    nrow = 1000,
                    ncol = nrow(.data)
  )
  
  probabilities <- .data[,"strength"] / sum(.data[,"strength"])
  # number of samples equals the duration of recall period divided by how long each sample takes (i.e., sampling rate)

  for (sim in 1:nrow(sim_RT)) {
    
    # The s vector tells you which item was sampled at each point in time
    s <- sample(x = .data[, "item"],
                size = duration,
                replace = TRUE,
                prob = probabilities)
    
    # This picks out the points in time where a new item was recalled
    RTs <- which(!duplicated(s))
    # This picks out the actual items output at those points in time, in the order they were output
    items <- s[RTs]
    # This stores the RTs for those items in this simulation
    sim_RT[sim, items] <- RTs
    
  }
  
  return(sim_RT)
}