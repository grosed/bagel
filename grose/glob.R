 mean <- function (x, ...) {
   UseMethod("mean", x)
 }

mean.numeric <- function(x, ...) sum(x) / length(x)
mean.data.frame <- function(x, ...) sapply(x, mean, ...)
mean.matrix <- function(x, ...) apply(x, 2, mean)mean.numeric <- function(x, ...) sum(x) / length(x)
