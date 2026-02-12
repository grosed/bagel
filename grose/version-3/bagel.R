
# type
setClass("bagel_kv_type", slots=list(prior = "function", H = "function", p0 = "numeric", p = "numeric", sigma = "numeric", particles = "list", t = "numeric"))
# constructor
bagel_kv <- function(prior,H,p0,p,sigma)
{
   return(new("bagel_kv_type", prior = prior, H = H, p0 = p0, p = p, sigma = sigma, particles = list(), t = 0L))
}


# update bagel object
setGeneric("update",function(object,y) standardGeneric("update"))
setMethod("update",c("bagel_kv_type","numeric"),
function(object,y)
{

   object@t <- object@t + 1L
   if(object@t == 1L)
   {
      # first data 
      tau <- 0L
      p <- particle_kv(object@prior,object@H,tau,object@p0,object@p,object@sigma)
      p <- theorem_2(p,object@t)
      p <- theorem_1(p,p,object@t)
      p <- theorem_4(p,object@t,y)
      p <- theorem_3(p,object@t,y)
      object@particles <- append(object@particles,p)
   }
   else
   {
      object@particles <- append(object@particles,theorem_2(object@particles[[1]],object@t)) # add the new particle
      object@particles <- Map(theorem_1 %><% list(object.1 = object@particles[[1]],t = object@t),object@particles)
      object@particles <- Map(theorem_4 %><% list(t = object@t,y = y),object@particles)
      object@particles <- Map(theorem_3 %><% list(t = object@t,y = y),object@particles)
      sumW <- Reduce("+",Map(function(.) .@weight,object@particles),0)
      object@particles <- Map(function(p) return(set_weight(p,p@weight/sumW)),object@particles)
   }
   return(object)
})


# update bagel object
setGeneric("weights",function(object) standardGeneric("weights"))
setMethod("weights",c("bagel_kv_type"),
function(object)
{
   weights <- unlist(Map(function(.) .@weight,object@particles))
   return(weights)
})












