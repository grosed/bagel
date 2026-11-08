library(bagelR)

result <- function(w, log_w0, threshold = Inf,
                   n = length(log_w0)) {
  bagel_result(y = rep(0, n), threshold = threshold,
               max_particles = length(w), w0ts = exp(log_w0),
               weights = w, posteriors = list(),
               ratios = as.list(rep(1, length(w))),
               taus = seq_along(w) - 1, t = length(log_w0) + 1,
               trace = list(), log_w0ts = log_w0)
}

full <- result(c(.2, .3, .5), log(c(1, .8, .2)))
stopifnot(identical(changepoint(full), list(location = 2L, detected = NA_integer_)))
no_alarm <- result(c(.8, .15, .05), log(c(1, .9, .8)), threshold = .95)
stopifnot(changepoint(no_alarm)$location == 1L,
          is.na(changepoint(no_alarm)$detected))

early <- result(c(.01, .09, .9), log(c(1, .8, .01)), threshold = .95, n = 5)
final <- result(c(.01, .09, .9), log(c(1, .8, .01)), threshold = .95)
stopifnot(identical(changepoint(early), list(location = 2L, detected = 3L)),
          identical(changepoint(final), changepoint(early)))

extreme <- result(c(0, 1), c(0, -1000), threshold = .95)
boundary <- result(c(.5, .5), log(c(1, .5)), threshold = .5)
no_early_stop <- result(c(0, 1), c(0, -1000), threshold = 1)
stopifnot(changepoint(extreme)$detected == 2L,
          is.na(changepoint(boundary)$detected),
          is.na(changepoint(no_early_stop)$detected))

single <- result(1, 0)
zero <- result(c(1, 0, 0), c(0, 0, 0))
stopifnot(is.na(changepoint(single)$location),
          is.na(changepoint(single)$detected),
          is.na(changepoint(zero)$location))
cat("Changepoint accessor regression tests passed.\n")
