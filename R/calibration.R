#' A Reference Class to generates different \code{calibrate.class} objects
#'
#' @description See the function \code{\link{calibrate}} which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link{calibrate}}... Other methods
#' should not be called as they are designed to be used during the calibration process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field md a \code{\link{model.class}} object generated by \code{\link{model}}
#' @field pr a \code{\link{prior.class}} object generated by \code{\link{prior}}
#' @field opt.estim list of estimation options (see \code{\link{calibrate}} for details)
#' @field opt.valid list of cross validation options (see \code{\link{calibrate}} for details)
#' @field logPost the log posterior
#' @field mcmc a \code{coda} variable of the generated chains
#' @field output list of several chains and acceptation ratios
#' @field onlyCV activates only the cross validation
#' @field errorCV output of the CV errors
#' @field n.cores number cores available
#' @field binf the lower bound of the parameters for the DOE
#' @field bsup the upper bound of the parameters for the DOE
#' @export
calibrate.class <- R6Class(classname = "calibrate.class",
                           public = list(
                             md         = NULL,
                             pr         = NULL,
                             opt.estim  = NULL,
                             opt.valid  = NULL,
                             logPost    = NULL,
                             mcmc       = NULL,
                             output     = NULL,
                             onlyCV     = NULL,
                             errorCV    = NULL,
                             ResultsCV  = NULL,
                             q05        = NULL,
                             q95        = NULL,
                             coverRate  = NULL,
                             n.cores    = NULL,
                             binf       = NULL,
                             bsup       = NULL,
                             initialize = function(md=NA,pr=NA,opt.estim=NA,opt.valid=NULL,onlyCV)
                             {
                               self$onlyCV <- onlyCV
                               self$md         <- md
                               self$pr         <- pr
                               self$opt.estim  <- opt.estim
                               self$opt.valid  <- opt.valid
                               self$logPost    <- private$logLikelihood(self$md$model)
                               Dim             <- length(self$pr)
                               if (Sys.info()[['sysname']]=="Windows")
                               {
                                 self$n.cores    <- 1
                               } else
                               {
                                 self$n.cores    <- 2
                               }
                               if (self$opt.estim$burnIn > self$opt.estim$Nmh)
                               {
                                 stop("The burnIn must be inferior to Nmh")
                               }
                               if (self$onlyCV==FALSE)
                               {
                                 if (opt.estim$Nchains==1)
                                 {
                                   self$output  <- self$calibration()
                                   self$mcmc    <- as.mcmc(self$output$out$THETA)
                                   chain        <- self$output$out$THETA[-c(1:self$opt.estim$burnIn),]
                                   qq           <- private$quantiles(chain)
                                   self$q05     <- qq$q05
                                   self$q95     <- qq$q95
                                 } else
                                 {
                                   cat("\nMultichain calibration is parallelized...\n\n")
                                   n            <- c(1:opt.estim$Nchains)
                                   if (Sys.info()[['sysname']]=="Windows")
                                   {
                                     self$output  <- lapply(n,self$calibration)
                                   } else
                                   {
                                     self$output  <- mclapply(n,self$calibration,mc.cores = self$n.cores)
                                   }
                                   self$mcmc    <- list()
                                   for (i in 1:opt.estim$Nchains)
                                   {
                                     self$mcmc[[i]] <- as.mcmc(self$output[[i]]$out$THETA)
                                   }
                                   self$output    <- self$output[[1]]
                                   chain          <- self$output$out$THETA[-c(1:self$opt.estim$burnIn),]
                                   qq             <- private$quantiles(chain)
                                   self$q05       <- qq$q05
                                   self$q95       <- qq$q95
                                   self$mcmc      <- as.mcmc.list(self$mcmc)
                                 }
                                 cat("\nEnd of the regular calibration\n\n")
                               }
                               if (is.null(self$opt.valid)==FALSE)
                               {
                                 cat(paste("\nThe cross validation is currently running on your ",
                                             self$n.cores," cores available....\n",sep=""))
                                 if (Sys.info()[[1]]=="Windows")
                                 {
                                   Results <- lapply(c(1:opt.valid$nCV),
                                                       self$CV)
                                 } else
                                 {
                                   Results <- mclapply(c(1:opt.valid$nCV),
                                                       self$CV,mc.cores = self$n.cores)
                                 }
                                 self$coverRate <- 0
                                 for (i in 1:self$opt.valid$nCV)
                                 {
                                   TempRes <- Results[[i]]
                                   inc <- TempRes$inc
                                   if (i==1){
                                     self$ResultsCV <- data.frame(Predicted=TempRes$Predict,Real=TempRes$Yval,Error=TempRes$err)
                                     if (self$onlyCV==FALSE)
                                     {
                                       if (TempRes$Predict < self$q95[inc] & TempRes$Predict > self$q05[inc])
                                       {
                                         self$coverRate <- self$coverRate + 1
                                       }
                                     }
                                     } else{
                                       self$ResultsCV <- rbind(self$ResultsCV,data.frame(Predicted=TempRes$Predict,
                                                                                         Real=TempRes$Yval,Error=TempRes$err))
                                       if (self$onlyCV==FALSE)
                                       {
                                         if (TempRes$Predict < self$q95[inc] & TempRes$Predict > self$q05[inc])
                                         {
                                           self$coverRate <- self$coverRate + 1
                                         }
                                       }
                                     }
                                 }
                                 self$coverRate <- self$coverRate/self$opt.valid$nCV
                                 self$errorCV <- sqrt(mean(self$ResultsCV$Error))
                               }
                             },
                             calibration = function(i=NA)
                             {
                               self$binf     <- private$boundaries()$binf
                               self$bsup     <- private$boundaries()$bsup
                               MetropolisCpp <- private$MCMC(self$md$model)
                               out           <- MetropolisCpp(self$opt.estim$Ngibbs,self$opt.estim$Nmh,
                                                              self$opt.estim$thetaInit,self$opt.estim$r,
                                                              self$opt.estim$sig,self$md$Yexp,self$binf,self$bsup,self$logPost,1)
                               MAP           <- private$MAPestimator(out)
                               return(list(out=out,MAP=MAP))
                             },
                             CV = function(i=NA)
                             {
                               Dim <- length(self$pr)
                               if (self$opt.valid$type.valid=="loo")
                               {
                                 inc <- sample(c(1:self$md$n),1)
                                 if (is.matrix(self$md$X))
                                 {
                                   dataCal <- self$md$X[-inc,]
                                   dataVal <- as.matrix(t(self$md$X[inc,]))
                                 } else
                                 {
                                   dataCal <- self$md$X[-inc]
                                   dataVal <- self$md$X[inc]
                                 }
                                 Ycal        <- self$md$Yexp[-inc]
                                 Yval        <- self$md$Yexp[inc]
                                 mdTempCal   <- model(code = self$md$code, X = dataCal, Yexp = Ycal,
                                                      model = self$md$model,
                                                      opt.gp = self$md$opt.gp,
                                                      opt.emul = self$md$opt.emul,
                                                      opt.sim = self$md$opt.sim)
                                 mdTempfit   <- self$calibrationCV(mdTempCal,Ycal)
                                 thetaMAP    <- mdTempfit$MAP
                                 if (mdTempCal$model=="model1" | mdTempCal$model=="model2")
                                 {
                                    Predict  <- mdTempCal$model.fun(thetaMAP[1:(Dim-1)],thetaMAP[Dim],dataVal)$y
                                    if (is.na(Predict)){Predict <- Yval}
                                 } else
                                 {
                                    Predict  <- mdTempCal$model.fun(thetaMAP[1:(Dim-3)],thetaMAP[(Dim-2):(Dim-1)],
                                                               thetaMAP[Dim],dataVal)$y
                                    if (is.na(Predict)){Predict <- Yval}
                                 }
                                 err <- sqrt((Predict-Yval)^2)
                                 res <- list(Predict=Predict,Yval=Yval,err=err,inc=inc)
                                 return(res)
                               }
                             },
                             calibrationCV = function(mdTemp,y)
                             {
                               binf          <- private$boundaries()$binf
                               bsup          <- private$boundaries()$bsup
                               MetropolisCpp <- private$MCMC(mdTemp$model)
                               out           <- MetropolisCpp(self$opt.estim$Ngibbs,self$opt.estim$Nmh,
                                                              self$opt.estim$thetaInit,self$opt.estim$r,
                                                              self$opt.estim$sig,self$md$Yexp,binf,bsup,self$logPost,0)
                               MAP           <- private$MAPestimator(out)
                               return(list(out=out,MAP=MAP))
                             }
                           ))


calibrate.class$set("private","MCMC",
                  function(model)
                  {
                    switch(model,
                           model1={return(MetropolisHastingsCpp)},
                           model2={return(MetropolisHastingsCpp)},
                           model3={return(MetropolisHastingsCppD)},
                           model4={return(MetropolisHastingsCppD)}
                    )
                  })


calibrate.class$set("private","boundaries",
                function()
                {
                  binf <- self$pr[[1]]$binf
                  bsup <- self$pr[[1]]$bsup
                  for (i in 2:length(self$pr))
                  {
                    binf <- c(binf,self$pr[[i]]$binf)
                    bsup <- c(bsup,self$pr[[i]]$bsup)
                  }
                  return(list(binf=binf,bsup=bsup))
                })


calibrate.class$set("private","quantiles",
                    function(chain)
                    {
                      Dist <- matrix(nr=nrow(chain),nc=length(self$md$Yexp))
                      dim   <- length(self$pr)
                      if (self$md$model == "model1" || self$md$model == "model2")
                      {
                        parFun <- function(i)
                        {
                          D  <- self$md$model.fun(chain[i,1:(dim-1)],chain[i,dim],self$md$X)$y
                          return(D)
                        }
                        if (Sys.info()[['sysname']]=="Windows")
                        {
                          res <- lapply(1:nrow(chain),parFun)
                        } else
                        {
                          res <- mclapply(1:nrow(chain),parFun,mc.cores = self$n.cores)
                        }
                        for (i in 1:nrow(chain))
                        {
                          Dist[i,] <- res[[i]]
                        }
                        qq <- apply(Dist,2,quantile,probs=c(0.05,0.95))
                        return(list(q05=qq[1,],q95=qq[2,]))
                      } else
                      {
                        parFun <- function(i)
                        {
                          D  <- self$md$model.fun(chain[i,1:(dim-3)],chain[i,(dim-2):(dim-1)],chain[i,dim],
                                                  self$md$X,CI=NULL)$y
                          return(D)
                        }
                        if (Sys.info()[['sysname']]=="Windows")
                        {
                          res <- lapply(1:nrow(chain),parFun)
                        } else
                        {
                          res <- mclapply(1:nrow(chain),parFun,mc.cores = self$n.cores)
                        }
                        for (i in 1:nrow(chain))
                        {
                          Dist[i,]  <- res[[i]]
                        }
                        qq <- apply(Dist,2,quantile,probs=c(0.05,0.95))
                        return(list(q05=qq[1,],q95=qq[2,]))
                      }
                    })


calibrate.class$set("private","logLikelihood",
                function(model)
                {
                  switch(model,
                         model1={return(private$logTest)},
                         model2={return(private$logTest)},
                         model3={return(private$logTestD)},
                         model4={return(private$logTestD)}
                  )
                }
)

calibrate.class$set("private","logTest",
                    function(theta,sig2)
                    {
                      if (length(self$pr) == 1)
                      {
                        return(log(self$md$likelihood(theta,sig2))+self$pr$prior(theta))
                      } else
                      {
                        s <- 0
                        for (i in 1:(length(theta)))
                        {
                          s <- s + self$pr[[i]]$prior(theta[i])
                        }
                        s <- s + self$pr[[(length(theta)+1)]]$prior(sig2)
                        return(self$md$likelihood(theta,sig2) + s)
                      }
                    })

calibrate.class$set("private","logTestD",
                function(theta,thetaD,sig2)
                {
                  s <- 0
                  for (i in 1:(length(theta)))
                  {
                    s <- s + self$pr[[i]]$prior(theta[i])
                  }
                  for (j in 1:(length(thetaD)))
                  {
                    s <- s + self$pr[[length(theta)+j]]$prior(thetaD[j])
                  }
                  s <- s + self$pr[[(length(theta)+1)]]$prior(sig2)
                  return(self$md$likelihood(theta,thetaD,sig2) + s)
                })

calibrate.class$set("private","MAPestimator",
                    function(out)
                      {
                        chain <- out$THETA[-c(1:self$opt.estim$burnIn),]
                        dens <- apply(chain,2,density)
                        map <- function(dens)
                        {
                          dens$x[which(dens$y==max(dens$y))]
                        }
                        return(unlist(lapply(dens,map)))
                    })



calibrate.class$set("public","plot",
                    function(x,graph="all",label=NULL,...)
                    {
                      if (missing(x)) stop("No x-axis selected, no graph is displayed",call. = FALSE)
                      p <- ncol(self$output$out$THETA)-1
                      p1 <- p2 <- p3 <- P <- L <- corrplot <- list()
                      inc = 1
                      for (i in 1:(p+1))
                      {
                        if (i == (p+1)) NameParam <- expression(sigma[err]^2)
                        else if (self$md$model %in% c("model3","model4") & i == p) NameParam <- expression(Psi[d])
                        else if (self$md$model %in% c("model3","model4") & i == (p-1)) NameParam <- expression(sigma[d]^2)
                        else if (!is.null(label)) NameParam <- label[i]
                        else NameParam <- substitute(theta[i])
                        p1[[i]] <- self$acf(i)+xlab("Lag")+ylab("ACF")
                        p2[[i]] <- self$mcmcChains(i) + xlab("iterations") + ylab(NameParam)
                        dplot2  <- data.frame(data=self$output$out$THETA[-c(1:self$opt.estim$burnIn),i],
                                              type="posterior")
                        p3[[i]] <- self$pr[[i]]$plot()+geom_density(data=dplot2,kernel="gaussian",adjust=3,alpha=0.1)+
                          xlab(NameParam)+theme(legend.text=element_text(size = 10),
                                                axis.text=element_text(size=20), axis.title = element_text(size=20))+
                          ylab("density")
                        if (i < (p+1))
                        {
                          for (j in 1:p)
                          {
                            if (j != i)
                            {
                              corrplot[[inc]] <- self$corrPlot(i,j)+ xlab(substitute(theta[i])) + ylab(substitute(theta[j]))
                              L[[j]] <- corrplot[[inc]]
                              inc <- inc+1
                            } else
                            {
                              L[[j]] <- p3[[i]]
                            }
                          }
                          arrangeGrob2 <- function(...) return(arrangeGrob(...,ncol = 1))
                          if (p != 1) P[[i]] <- do.call(arrangeGrob2,L)
                        }
                      }
                      output <- self$outputPlot(select.X=x)+ xlab("x") +ylab("y")
                      grid.arrange2 <- function(...) return(arrangeGrob(...,nrow=1))
                      chainPlotLayout <- function()
                      {
                        if (self$md$model %in% c("model1","model2")) N <- p+1 else N <- p+3
                        t  <- N%/%3
                        tr <- N%%3
                        inc <- 0
                        if (t == 0)
                        {
                          grid.arrange(do.call(grid.arrange2,p1),do.call(grid.arrange2,p2),do.call(grid.arrange2,p3))
                        } else
                        {
                          for (j in 1:t)
                          {
                            inc <- inc + 1
                            if (inc > N) break
                            grid.arrange(do.call(grid.arrange2,p1[(3*(j-1)+1):(3*j)]),
                                         do.call(grid.arrange2,p2[(3*(j-1)+1):(3*j)]),
                                         do.call(grid.arrange2,p3[(3*(j-1)+1):(3*j)]))
                          }
                        }
                      }
                      corrPlotLayout <- function()
                      {
                        if (p == 1) {warning("No correlation plot available",call. = FALSE)}
                        else
                        {
                           plot(do.call(grid.arrange2,P))
                        }
                      }
                      if (!is.null(graph))
                      {
                        if ("chains" %in% graph)
                        {
                          chainPlotLayout()
                        } else if ("corr" %in% graph)
                        {
                          corrPlotLayout()
                        } else if ("result" %in% graph)
                        {
                          print(output)
                        } else if (graph == "all")
                        {
                          chainPlotLayout()
                          corrPlotLayout()
                          print(output)
                        } else
                        {
                          warning("The graph value is not correct all graphs are displayed",call. = FALSE)
                          chainPlotLayout()
                          corrPlotLayout()
                          print(output)
                        }
                      } else {
                        #warning('Nothing is plot since the graph option is NULL',call. = FALSE)
                      }
                      invisible(list(ACF=p1,MCMC=p2,dens=p3,corrplot=corrplot,out=output))
                    })


calibrate.class$set("public","acf",
                    function(i)
                    {
                      bacf   <- acf(self$output$out$THETA[-c(1:self$opt.estim$burnIn),i], plot = FALSE)
                      bacfdf <- with(bacf, data.frame(lag, acf))
                      p      <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf))+
                        geom_hline(aes(yintercept = 0))+
                        geom_segment(mapping = aes(xend = lag, yend = 0))+
                        xlab("")+ylab("")+theme_light()
                      return(p)
                      })


calibrate.class$set("public","corrPlot",
                    function(i,j)
                    {
                      chain1 <- self$output$out$THETA[-c(1:self$opt.estim$burnIn),i]
                      chain2 <- self$output$out$THETA[-c(1:self$opt.estim$burnIn),j]
                      df     <- data.frame(x=chain1,y=chain2)
                      p      <- ggplot(data = df, mapping = aes(x=x, y=y))+
                        geom_jitter() + xlab(expression(theta_1))+
                        ylab(expression(theta_1))+theme_light()
                      return(p)
                    })

calibrate.class$set("public","mcmcChains",
                    function(i)
                    {
                      n <- length(self$output$out$THETA[-c(1:self$opt.estim$burnIn),i])
                      resgg <- data.frame(inc=c(1:n),data=self$output$out$THETA[-c(1:self$opt.estim$burnIn),i])
                      p   <- ggplot(data=resgg, aes(x=inc,y=data))+geom_line()+ylab("")+
                        xlab("")+theme_light()
                      return(p)
                    }
                    )


calibrate.class$set("public","outputPlot",
                    function(select.X=NULL)
                    {
                      if (is.null(select.X)==TRUE)
                        {
                          if (is.null(dim(self$md$X))==TRUE)
                          {
                            X <- self$md$X
                          }else
                          {
                            stop('The dimension of X is higher than 1, the plot cannot be provided for a dimension >1')
                          }
                        } else {X <- select.X}
                      m <- self$output$out$THETA[-c(1:self$opt.estim$burnIn),]
                      Dist <- Dist2 <- matrix(nr=nrow(m),nc=length(self$md$Yexp))
                      dim   <- length(self$pr)
                      if (self$md$model == "model1" || self$md$model == "model2")
                      {
                        parFun <- function(i)
                        {
                          D  <- self$md$model.fun(m[i,1:(dim-1)],m[i,dim])$y
                          return(D)
                        }
                        if (Sys.info()[['sysname']]=="Windows")
                        {
                          res <- lapply(1:nrow(m),parFun)
                        } else
                        {
                          res <- mclapply(1:nrow(m),parFun,mc.cores = self$n.cores)
                        }
                        for (i in 1:nrow(m))
                        {
                          Dist[i,] <- res[[i]]
                        }
                      } else
                      {
                        print('The computational time might be long to get the output plot')
                        parFun <- function(i)
                        {
                          D  <- self$md$model.fun(m[i,1:(dim-3)],m[i,(dim-2):(dim-1)],m[i,dim])$y
                          D2 <- self$md$model.fun(self$output$MAP[1:(dim-3)],m[i,(dim-2):(dim-1)],
                                                   self$output$MAP[dim])$y
                          return(list(D=D,D2=D2))
                        }
                        if (Sys.info()[['sysname']]=="Windows")
                        {
                          res <- lapply(1:nrow(m),parFun)
                        } else
                        {
                          res <- mclapply(1:nrow(m),parFun,mc.cores = self$n.cores)
                        }
                        for (i in 1:nrow(m))
                        {
                          Dist[i,]  <- res[[i]]$D
                          Dist2[i,] <- res[[i]]$D2
                        }
                        qqd <- apply(Dist2,2,quantile,probs=c(0.05,0.95))
                        ggdata2 <- data.frame(y=self$md$Yexp,x=X,upper=qqd[2,],lower=qqd[1,],type="experiments",
                                            fill="95% credibility interval for the discrepancy")
                      }
                      qq <- apply(Dist,2,quantile,probs=c(0.05,0.95))
                      if (self$md$model=="model1"||self$md$model=="model2")
                      {
                        MAP <- self$output$MAP
                        p <- length(MAP)-1
                        Ys <- self$md$model.fun(MAP[1:p],MAP[p+1])$y
                      }else
                      {
                        MAP <- self$output$MAP
                        p <- length(MAP)-3
                        Ys <- self$md$model.fun(MAP[1:p],MAP[(p+1):(p+2)],MAP[p+3])$y
                      }
                      ggdata <- data.frame(y=Ys,x=X,upper=qq[2,],lower=qq[1,],type="calibrated",
                                           fill="95% credibility interval a posteriori")
                      if (self$md$model == "model1" || self$md$model == "model2")
                      {
                        ggdata2 <- data.frame(y=self$md$Yexp,x=X,upper=qq[2,],lower=qq[1,],type="experiments",
                                              fill="95% credibility interval a posteriori")
                      }
                      ggdata <- rbind(ggdata,ggdata2)
                      p <- ggplot(ggdata) + geom_line(aes(x=x,y=y,color=type))+
                        geom_ribbon(aes(ymin=lower, ymax=upper, x=x,fill=fill), alpha = 0.3) +
                        scale_fill_manual("",values=c("blue4","black")) +
                        theme_light() + xlab("") + ylab("") +
                        theme(legend.title = element_blank(),
                              legend.position = c(0.2,0.8),
                              legend.background = element_rect(
                                linetype="solid", colour ="grey"))
                     return(p)
                    })


calibrate.class$set("public","print",
                    function()
                    {
                      cat("Call:\n\nWith the function:\n")
                      print(self$md$code)
                      cat("\nSelected model :",self$md$model,"\n")
                      if (self$onlyCV==TRUE)
                      {
                        cat("\nCross validation:\n Method:",self$opt.valid$type.valid,"\n")
                        print(head(self$ResultsCV))
                        cat("\nRMSE: ")
                        print(self$errorCV)
                      } else
                      {
                        cat("\nAcceptation rate of the Metropolis within Gibbs algorithm:\n")
                        print(paste(round(self$output$out$AcceptationRatioWg/self$opt.estim$Ngibbs*100,2),
                                    "%",sep = ""))
                        cat("\nAcceptation rate of the Metropolis Hastings algorithm:\n")
                        print(paste(self$output$out$AcceptationRatio/self$opt.estim$Nmh*100,"%",sep = ""))
                        cat("\nMaximum a posteriori:\n")
                        print(apply(self$output$out$THETA[-c(1:self$opt.estim$burnIn),],2,max))
                        cat("\nMean a posteriori:\n")
                        print(apply(self$output$out$THETA[-c(1:self$opt.estim$burnIn),],2,mean))
                        cat("\n")
                        if (is.null(self$opt.valid)==FALSE)
                        {
                          cat("\nCross validation:\n Method:",self$opt.valid$type.valid,"\n")
                          print(head(self$ResultsCV))
                          cat("\nRMSE: ")
                          print(self$errorCV)
                          if (self$onlyCV==FALSE)
                          {
                            cat("\nCover rate: \n")
                            print(paste(round(self$coverRate*100,2),"%",sep = ""))
                          } else
                          {
                            cat("\nTo activate the coverage rate, please desable the otion onlyCV")
                          }

                        }
                      }
                    }
)

