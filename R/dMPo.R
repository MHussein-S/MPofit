dMPo<-function (x,par,distr,log = FALSE)
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
  g<-do.call(ddistname, c(list(x), as.list(distpar)))
  G<-do.call(pdistname, c(list(x), as.list(distpar)))
  d<-a^(G-1)*g*(1+G*log(a))
  if (log)
    d <- log(d)
  return(d)

}
