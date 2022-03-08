#cdf
pMPo<-function (q,par,distr, lower.tail = TRUE, log.p = FALSE )
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
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
  G<-do.call(pdistname, c(list(q), as.list(distpar)))
  p<-a^(G-1)*G
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  return(p)
}
