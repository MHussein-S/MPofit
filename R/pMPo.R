#cdf
pMPo<-function (q,par,distr, lower.tail = TRUE, log.p = FALSE )
{
  if(is.null(names(par)))
    stop(" 'par' must be a named list.")
  ddistname <- paste("d", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  argdistname <- names(formals(ddistname))
  if (!exists(ddistname, mode = "function"))
    stop("The ", ddistname, " function must be defined")
  if (!is.list(par))
    par<-as.list(par)
  #exclude a from par
  apar<-match("a",names(par))
  if(!is.na(apar))
  {
    a<-par$a
    par<-par[-apar] #exclude a parameter from par
  }
  else
    stop(" 'a' parameter not defined")
  m <- match(names(par), argdistname)
  if (any(is.na(m)))
  {
    stop(" 'par' must specify names which are arguments to ", distr)
  }
  G<-do.call(pdistname, c(list(q), as.list(par)))
  p<-a^(G-1)*G
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  return(p)
}
