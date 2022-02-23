#Quantile function of the MPo distribution
qMPo<-function(p,par,distr,lower=0,upper,lower.tail=TRUE,log.p=FALSE)
{
  if(is.null(names(par))|any(names(par)==""))
    stop("'par' must be a named numeric vector of the form 'c(name=val,name=val,...)'")
  if (!is.character(distr))
    stop("distr must be a character string naming the baseline distribution")

  setpar<-function(ddistname,allpar)
  {
    if(!is.list(allpar))
      allpar<- as.list(allpar)
    if (!exists(ddistname, mode="function"))
      stop(paste("The ", ddistname, " distribution is not defined"))
    args <- names(formals(ddistname))
    a<-allpar$a
    extrpar<-list(a=a)
    distparn<-setdiff(names(allpar),c("a"))
    distpar<-allpar[distparn]
    m <- match(distparn,args)
    if (any(is.na(m)))
      stop("you specifies names of parameters which are not valid for ",ddistname)
    return(list(extrpar=extrpar,distpar=distpar))
  }
  ddistname <- paste("d", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  parset<-setpar(ddistname,allpar=par)
  a<-parset$extrpar$a
  distpar<-parset$distpar
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


