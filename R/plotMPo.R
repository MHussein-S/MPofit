#Empirical and theoretical: dens, CDF
#and P-P plot
plotMPo<-function (data, distr, para, histo = TRUE, breaks = "default", demp = TRUE)
{
  if (missing(data) || !is.vector(data, mode = "numeric"))
    stop("data must be a numeric vector")
  if ((missing(distr) & !missing(para)) || (!missing(distr) &  missing(para)))
    stop("distr and para -must defined")
  if (!histo & !demp)
    stop("one the arguments histo and demp must be put to TRUE")
  s <- sort(data)
  obsp<-ppoints(s) #observed probabilities
  if (length(s) != length(obsp))
    stop("problem when computing probabilities.")
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
  n <- length(data)
  if (missing(distr))
  {
    par(mfrow = c(1, 2))
    if (histo)
    {
      if (demp)
      {
        if (breaks == "default")
          h <- hist(data, freq = FALSE, xlab = "Data", main = "Empirical density")
        else
          h <- hist(data, freq = FALSE, xlab = "Data", main = "Empirical density", breaks = breaks)
        lines(density(data)$x, density(data)$y, lty = 1, col	= "red")
      }
      else
      {
        if (breaks == "default")
          h <- hist(data, freq = FALSE, xlab = "Data", main ="Histogram")
        else
          h <-hist(data, freq = FALSE, xlab = "Data", main ="Histogram", 	breaks	= breaks)
      }
    }
    else
    {
      h <- hist(data, freq = FALSE, plot = FALSE)

      plot(density(data)$x, density(data)$y, lty = 1, col = "red", type = "l", xlab = "Data", main = paste("Empirical density"), ylab = "Density")

    }
    plot(s, obsp, main = paste("Cumulative distribution"), xlab = "Data", xlim = c(min(h$breaks), max(h$breaks)), ylab = "CDF")
  }
  else #distr not missing
  {
    if (is.null(names(para)))
      stop("'para' must be a named list or vector")
    ddistname <- paste("d", distr, sep = "")
    pdistname <- paste("p", distr, sep = "")
    parset<-setpar(ddistname,allpar=as.list(para))
    a<-parset$extrpar$a
    distpar<-parset$distpar
    par(mfrow = c(2, 2))
    if (breaks == "default")
      h <- hist(data, plot = FALSE)
    else
    h <- hist(data, breaks = breaks, plot = FALSE)
    xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)

    G<-do.call(pdistname, c(list(xhist), as.list(distpar)))
    g<-do.call(ddistname, c(list(xhist), as.list(distpar)))
    F<-a^(G-1)*G
    f<-a^(G-1)*g*(1+G*log(a))
    yhist <- f
    if (length(yhist) != length(xhist))
      stop("problem when computing densities.")
    ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(max(h$density),max(density(data)$y)))
    if (histo)
    {
      hist(data, freq = FALSE, xlab = "Data", ylim = c(0, ymax), breaks = h$breaks, main = paste("Empirical and theoretical dens."),ylab= "Density", xlim = c(min(h$breaks), max(h$breaks)))
      lines(xhist,yhist , lty = 2, col = "black", type = "l", xlab = "Data", main = paste("Empirical and theoretical dens."), ylab= "Density", 	xlim = c(min(h$breaks), max(h$breaks)))
      if (demp)
      {
       lines(stats::density(data)$x, stats::density(data)$y, lty = 1, col = "red")
      legend("topright", bty = "n", lty = c(2,1), col = c("black", "red"), legend = c("theoretical","empirical"), bg = "white", cex = 0.7)
      }
    }
    else
    {
      plot(xhist,yhist , lty = 2, col = "black", type = "l", xlab = "Data", main = paste("Empirical and theoretical dens."), ylab= "Density", 	xlim = c(min(h$breaks), max(h$breaks)))
      if (demp)
      {
        lines(stats::density(data)$x, stats::density(data)$y, lty = 1, col = "red",xlim = c(min(h$breaks), max(h$breaks)))
        legend("topright", bty = "n", lty = c(2,1), col = c("black", "red"), legend = c("theoretical","empirical"), bg = "white", cex = 0.7)
      }
    }
    # Plot of the cumulative probability distributions
    xmin <- min(h$breaks)
    xmax <- max(h$breaks)
    plot(s, obsp, main = paste("Empirical and theoretical CDFs"),  xlab = "Data", 	ylab = "CDF", xlim = c(xmin, xmax)) #Empirical CDF

    sfin <- seq(xmin, xmax, by = (xmax - xmin)/100)
    #G_fin<-do.call(pdistname, c(list(sfin), as.list(distpar)))
    #F_fin<-a^(G_fin-1)*G_fin  #Theoretical CDF
    F_fin<-pMPo(sfin,par = para,distr = distr)  #Theoretical CDF
    lines(sfin, F_fin, lty = 1, col = "red")

    #P-P plot
  #  G_s<-do.call(pdistname, c(list(s), as.list(distpar)))
 #   F_s<-a^(G_s-1)*G_s

    F_s<-pMPo(s,par = para,distr = distr)  #Theoretical CDF
    if (length(F_s) != length(obsp))
      stop("problem when computing probabilities.")
    plot(F_s, obsp, main = "P-P plot", xlab = "Theoretical probabilities", ylab = "Empirical probabilities")
    abline(0, 1, col="red")
    #Q-Q plot
    theoQ<-qMPo(obsp,par=para,distr=distr,lower=0,upper=1000) #Theoretical quantiles
    if (length(theoQ) != length(obsp))
      stop("problem when computing quantities.")
    plot(theoQ, s, main = "Q-Q plot", xlab = "Theoretical quantiles", ylab = "Empirical quantiles")

    abline(0, 1, col="red")
  }
}


