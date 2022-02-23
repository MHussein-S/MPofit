# r-th ordinary moment of the new distribution for any base distribution
mMPo<-function(r,par,distr)
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " distribution is not defined"))
  apar<-match("a",names(par))
  if(is.na(apar))
    stop(" 'a' parameter is not defined")
  args <- names(formals(ddistname))
  a<-par$a
  distparn<-setdiff(names(par),c("a"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specifies names of parameters which are not valid for ",ddistname)
  if (a < exp(-1))
    stop("MPo distribution not defined for a < 1/e")

  ########################## integrand ##########################
  integrand <- function(y,a,qdistname,distpar,r)
  {
    L1<-do.call(qdistname, c(list(y), as.list(distpar)))
    L1<-a^y*(1+y*log(a))*L1^r
    L1<-L1/a
  }
  #-------------------------------------------
  result<-integrate(integrand, lower = 0, upper =1,a=a,qdistname=qdistname ,distpar=distpar,r=r)
  return(result$value)
}
#-------------------------------------------
