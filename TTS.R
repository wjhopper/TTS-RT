TTS <- function(x, ...) {
  UseMethod("TTS", x)
}

TTS.numeric <- function(strength_dist_properties = c(.5, .15),
                        list_length = 10,
                        threshold = 0,
                        duration = 30,
                        rate = 1) {

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
  
  if (length(rate) != 1 && length(rate) != list_length) {
    stop("'rate' argument must be a scalar or same length as the value of 'list_length' argument")
  }
  
  if (length(rate) == 1) {
    threshold <- rep(rate, list_length)
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
                       rate = 1) {
  
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
    # The rate argument scales the duration of the recall period
    # e.g., rate of 2 samples per second with a 30 second duration is the
    # same as a rate of 1 sample per second with a 60 second duration)
    s <- sample(x = .data[, "item"],
                size = duration*rate,
                replace = TRUE,
                prob = probabilities)
    
    item_index <- match(s, .data[,"item"])
    s_thresholds <- .data[item_index, "threshold"]
    s_strengths <- .data[item_index, "strength"]
    s[s_strengths <= s_thresholds] <- NA
    
    # This picks out the points in time where a new item was recalled
    RTs <- which(!duplicated(s, incomparables = NA) & !(is.na(s)))
    # This picks out the actual items output at those points in time, in the order they were output
    items <- s[RTs]
    # This stores the RTs for those items in this simulation
    sim_RT[sim, items] <- RTs
    
  }
  
  sim_RT <- sim_RT/rate # Rescale RT's back into seconds instead of "sampling-speed" units.
  
  x <- structure(list(RT = sim_RT,
                      memories = .data,
                      duration = duration,
                      rate = rate),
                 class = "TTS")
  return(x)
}

summary.TTS <- function(TTS_obj, fit = TRUE) {

  mean_RT_by_item <- colMeans(TTS_obj$RT, na.rm = TRUE)
  mean_acc_by_item <- apply(TTS_obj$RT, 2, function(x) mean(!is.na(x)))
  
  mean_RT <- mean(TTS_obj$RT, na.rm = TRUE)
  mean_acc <- mean(mean_acc_by_item)
  
  results_table <- rbind(c(mean_acc_by_item, mean_acc),
                         c(mean_RT_by_item, mean_RT)
                         )
  colnames(results_table) <- c(paste("Item", as.character(1:(ncol(results_table)-1))), "Avg.")
  rownames(results_table) <- c("Acc.", "R.T.")
  
  outputs_hist <- hist(TTS_obj$RT, breaks = seq(0,TTS_obj$duration), plot = FALSE)
  proportion_per_timepoint <- outputs_hist$density
  cum_proportion_per_timepoint <- cumsum(proportion_per_timepoint)*mean_acc

  x <- structure(list(means = results_table,
                      proportion_per_timepoint = proportion_per_timepoint,
                      cum_proportion_per_timepoint = cum_proportion_per_timepoint),
                 class = "TTS_summary"
  )

  # Fitting exponential distribution rate parameter to simulated data
  if (fit) {

    x$fitted_rate <- mle2(y ~ dexp(rate = lambda),
                          start = list(lambda = 1 / ncol(TTS_obj$RT)),
                          data = data.frame(y = TTS_obj$RT[!is.na(TTS_obj$RT)])
                          )
  }

  return(x)
}

cum_exp <- function(x, rate, asymptote) {
  asymptote * (1 - exp(-x*rate))
}

plot.TTS_summary <- function(TTS_obj){

  par(mfrow=c(1,2))
  
  # Latency Distibution, with best fitting exponential if
  plot(x = seq_along(TTS_obj$proportion_per_timepoint),
       y = TTS_obj$proportion_per_timepoint,
       xlab = "Recall period (s)", ylab = "Relative Frequency",
       main = "Recall Latencies")
  
  if ("fitted_rate" %in% names(TTS_obj)) {

    lines(x = c(0:length(TTS_obj$proportion_per_timepoint)),
          y = dexp(0:length(TTS_obj$proportion_per_timepoint),
                   rate = coef(TTS_obj$fitted_rate))
          )
    text(x = params$duration/2, y = max(TTS_obj$proportion_per_timepoint)*.8,
         labels = substitute(lambda == y, list(y=round(coef(TTS_obj$fitted_rate), 3)))
         )
  }
  
  # Cummulative Latency Distribution
  plot(x = seq_along(TTS_obj$cum_proportion_per_timepoint),
       y = TTS_obj$cum_proportion_per_timepoint,
       xlab = "Recall period (s)", ylab = "Proportion Recalled",
       ylim = c(0, 1), main = "Cummulative Latency"
       )
  curve(cum_exp(x,
                rate = coef(TTS_obj$fitted_rate),
                asymptote = TTS_obj$means[1, ncol(TTS_obj$means)]
                ),
        add = TRUE)
  
  par(mfrow=c(1,1))
}

coef.TTS_summary <- function(TTS_obj, ...) {
  return(coef(TTS_obj$fitted_rate, ...))
}

pander.TTS_summary <- function(TTS_obj, ...) {
  pander(TTS_obj$means, caption = "Average Performance", ...)
}