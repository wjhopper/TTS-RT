TTS <- function(x, ...) {
  UseMethod("TTS", x)
}

TTS.numeric <- function(strength_dist_properties = c(.5, .15),
                        list_length = 10,
                        threshold = 0,
                        duration = 30,
                        constant_rate = TRUE) {

  if (length(strength_dist_properties) != 2) {
    stop("'strength_dist_properties' argument must be a numeric vector with 2 elements")
  }
  
  if (length(list_length) != 1 || !is.numeric(list_length)) {
    stop("'list_length' argument must be a numeric scalar")
  }
  
  if (length(list_length) != 1 || !is.numeric(list_length)) {
    stop("'duration' argument must be a numeric scalar")
  }
  
  if (length(threshold) != 1 && length(threshold) != list_length) {
    stop("'threshold' argument must be a scalar or same length as the value of 'list_length' argument")
  }
  
  if (length(threshold) == 1) {
    threshold <- rep(threshold, list_length)
  }
  
  beta_parameters <- whoppeR::betaParams(mean = strength_dist_properties[1],
                                         sd = strength_dist_properties[2]
                                         )
  memory_strengths <- rbeta(list_length,
                            shape1 = beta_parameters$a,
                            shape2 = beta_parameters$b
                            )
  memory_set <- matrix(c(1:list_length, memory_strengths, threshold),
                       nrow = list_length,
                       dimnames = list(NULL, c("item", "strength", "threshold"))
                       )
  
  x <- TTS.matrix(.data = memory_set, duration = duration, constant_rate = constant_rate)
  return(x)
}

TTS.data.frame <- function(.data, ...) {
  
  is_column_numeric <- vapply(.data, is.numeric, FUN.VALUE = logical(1))
  if (!all(is_column_numeric)) {
    stop("All variables must be numeric")
  }

  x <- TTS.matrix(as.matrix(.data), ...)
  return(x)
}

TTS.matrix <- function(.data,
                       duration = 30,
                       constant_rate = TRUE) {
  
  if (!is.numeric(.data)) {
    stop("input matrix must be numeric")
  }
  
  if ( !all(c("item","strength","threshold") %in% colnames(.data)) ) {
    stop("input matrix must contain named columns 'item', 'strength', and 'threshold'")
  }
  
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