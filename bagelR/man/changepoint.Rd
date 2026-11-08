\name{changepoint}
\alias{changepoint}
\alias{changepoint,bagel_result_type-method}
\title{MAP changepoint location and alarm time}
\usage{
changepoint(object)
}
\arguments{
  \item{object}{A fitted \code{bagel_result_type} object.}
}
\description{
Returns the MAP changepoint location conditional on a change, together with
the first observation at which the stopping threshold was exceeded.
}
\value{
A list with \code{location} and \code{detected}. The location is available
even when no alarm occurred, including runs with \code{threshold = Inf}.
It is \code{NA_integer_} if there are no positive-weight change candidates.
Ties are resolved by the first maximum in the recovered location weights.
The alarm time is \code{NA_integer_} if the threshold was never exceeded.
An alarm at the final observation is included.
}
\details{
The alarm occurs when the posterior probability of a change exceeds
\code{threshold}, which is specified on the ordinary probability scale.
For example, \code{threshold = 0.95} means a change probability above 95\%.
A threshold of 1 or greater, including \code{Inf}, disables early stopping.
The comparison uses retained log probabilities internally for numerical
stability. The location uses all processed observations; it does not imply
that an alarm was raised.

This accessor always returns a list. Code that previously tested whether
the entire return value was \code{NA} should instead test
\code{is.na(changepoint(object)$detected)}.
}
