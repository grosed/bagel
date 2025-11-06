

kalman_filter <- function(P)  structure(P,class = kalman_filter)

update <- function(kf,a,b)
{
	UseMethod("update",a,b)
}

update.kalman_filter <- function(kf,a,b)
{
	print(a)
	print(b)
}


kf <- kalman_filter(5)

update(kf,3,4)

