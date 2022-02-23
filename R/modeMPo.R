# Mode of the new family by direct maximization of log(f(x))
modeMPo<-function(start=0,parm,distr,lower=0,upper)
  {
  if(is.null(names(parm))|any(names(parm)==""))
    stop("'parm' must be a named numeric vector of the form 'c(name=val,name=val,...)'")
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

  objfn <- function(param,x,distr)
    {
    ddistname <- paste("d",distr,sep="")
    pdistname <- paste("p",distr,sep="")
    parset<-setpar(ddistname,allpar=as.list(param))
    a<-parset$extrpar$a
    distpar<-parset$distpar
    G<-do.call(pdistname,c(list(x),as.list(distpar)))
    g<-do.call(ddistname,c(list(x),as.list(distpar)))
    Logg <- (G-1)*log(a)+log(g)+log(1+G*log(a))
    return(-Logg)
    }
  opt <- optim(par = start,fn=objfn, param = parm,distr=distr,method = "Brent",lower = lower,upper = upper)
  opt$par
  }
