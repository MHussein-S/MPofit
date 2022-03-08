#Quantile function of the MPo distribution
qMPo<-function(p,par,distr,lower=0,upper,lower.tail=TRUE,log.p=FALSE)
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")

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
  if (log.p==TRUE)
    p<-exp(p)
  if (lower.tail==FALSE)
    p<-1-p
  qfun<-function (x,p,a,distpar,pdistname)
  {
    G<-do.call(pdistname, c(list(x), as.list(distpar)))
    F<-a^(G-1)*G
    F-p
    }
  #--
  res<-c()
  for (pi in p)
    {
    uni <- uniroot(f=qfun, p=pi, a=a,distpar=distpar, pdistname=pdistname,f.lower = -pi,f.upper = pi,lower = lower,upper = upper)$root
    res<-c(res,uni)
    }
  return(res)
  }


