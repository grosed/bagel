# constructor
kalman_filter <- function(P,x,e = NA)
{
	kf <- list("P" = function() return(P),
		   "x" = function() return(x),
		   "e" = function() return(e))
        class(kf) <- "kalman_filter"
	return(kf)
}


kalman_filter_predict <- function(kf,F,B,Q,u)
{
	UseMethod("kalman_filter_predict")
}

kalman_filter_update <- function(kf,F,B,Q,H,R,u,z)
{
	UseMethod("kalman_filter_update")
}

kalman_filter_predict.kalman_filter <- function(kf,F,B,Q,u)
{
	P <- kf$P()
	x <- kf$x()
	return(list("x" = F%*%x + B%*%u,"P" = F%*%P%*%t(F) + Q))
}

kalman_filter_update.kalman_filter <- function(kf,F,B,Q,H,R,u,z)
{
	prediction <- kalman_filter_predict(kf,F,B,Q,u)
	x <- prediction$x
	P <- prediction$P
	y <- z - H%*%x
	S <- H%*%P%*%t(H) + R
	K <- P%*%H%*%solve(S)
	x <- x + K%*%y
	P <- P - K%*%H%*%P
	e <- z - H%*%x
	return(kalman_filter(P,x,e))
}



