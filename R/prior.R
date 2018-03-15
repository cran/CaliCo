#' A Reference Class to generates differents \code{\link{prior.class}} objects
#'
#'
#' @description See the function \code{\link{prior}} which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation \code{\link{prior}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field type.prior of the selected prior
#' @field opt.prior the charasteristics of the selected prior
#' @field log if we want the log result or not
#' @export
prior.class <- R6Class(classname = "prior.class",
                           public = list(
                             type.prior  = NULL,
                             opt.prior   = NULL,
                             log         = NULL,
                             initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                             {
                               self$type.prior  <- type.prior
                               self$opt.prior   <- opt.prior
                               self$log         <- log
                               private$checkPrior()
                             }
                           ))

prior.class$set("private","checkPrior",
                function()
                {
                  if (self$type.prior != "gaussian" && self$type.prior != "gamma" &&
                      self$type.prior != "invGamma" && self$type.prior != "unif")
                  {
                    stop('Please select a correct prior')
                  }
                })



gaussian.class <- R6Class(classname = "gaussian.class",
                              inherit = prior.class,
                              public = list(
                                mean  = NULL,
                                var   = NULL,
                                sd    = NULL,
                                binf  = NULL,
                                bsup  = NULL,
                                initialize = function(type.prior=NA,opt.prior=NA,log=TRUE)
                                {
                                  super$initialize(type.prior,opt.prior,log)
                                  self$mean  <- unlist(self$opt.prior)[1]
                                  self$var   <- unlist(self$opt.prior)[2]
                                  self$sd    <- sqrt(self$var)
                                  self$binf  <- self$mean - 4*sqrt(self$var)
                                  self$bsup  <- self$mean + 4*sqrt(self$var)
                                },
                                prior = function(x=seq(self$binf,self$bsup,length.out = 100))
                                  {
                                  if (self$log == FALSE)
                                  {
                                    return(1/(sqrt(2*pi*self$sd))*
                                               exp(-1/(2*self$sd)*
                                                     (x-self$mean)^2))
                                    } else
                                  {
                                    return(log(1/(sqrt(2*pi*self$sd))*
                                                   exp(-1/(2*self$sd)*
                                                         (x-self$mean)^2)))
                                  }
                                },
                                plot = function()
                                {
                                  dplot <- data.frame(data=rnorm(n=1000,self$mean,sqrt(self$var)),type="prior")
                                  p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                                    geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                                    theme_light()+xlab("")+ylab("")+
                                    theme(legend.position=c(0.86,0.86),
                                          legend.text=element_text(face="bold",size = '20'),
                                          legend.title=element_blank(),
                                          legend.key=element_rect(colour=NA),
                                          axis.text=element_text(size=20))+
                                    geom_hline(aes(yintercept = 0))
                                  return(p)
                                }
                              ))


unif.class <- R6Class(classname = "unif.class",
                              inherit = prior.class,
                              public = list(
                                binf    = NULL,
                                bsup    = NULL,
                                y       = NULL,
                                initialize = function(type.prior=NA,opt.prior=NA,log=TRUE)
                                {
                                  super$initialize(type.prior,opt.prior,log)
                                  self$binf <- unlist(self$opt.prior)[1]
                                  self$bsup <- unlist(self$opt.prior)[2]
                                },
                                prior = function(x=seq(self$binf,self$bsup,length.out = 100))
                                {
                                  self$y <- matrix(nr=length(x),1)
                                  for (i in 1:length(x))
                                  {
                                    if (x[i]>=self$binf & x[i]<=self$bsup)
                                    {
                                      self$y[i] <- 1/(self$bsup-self$binf)
                                    } else
                                    {
                                      self$y[i] <- 0
                                    }
                                  }
                                  if (self$log == FALSE)
                                  {
                                    return(self$y)
                                  } else
                                  {
                                    return(log(self$y))
                                  }
                                },
                                plot = function()
                                {
                                  xvals <- data.frame(data=c(self$binf,self$bsup),type="prior")
                                  ggplot(data.frame(xvals), aes(x = data,color=type)) +
                                    stat_function(fun = dunif, color="red")+
                                    stat_function(fun = dunif, args = list(min =self$binf,max =self$bsup)
                                                  , geom = "area",fill = "red", alpha = 0.1) +
                                    theme_light()+xlab("")+ylab("")+
                                    theme(legend.position=c(0.86,0.86),
                                          legend.text=element_text(face="bold",size = '20'),
                                          legend.title=element_blank(),
                                          legend.key=element_rect(colour=NA),
                                          axis.text=element_text(size=20))+
                                    geom_hline(aes(yintercept = 0))
                                }
                              ))


gamma.class <- R6Class(classname = "gamma.class",
                          inherit = prior.class,
                          public = list(
                            shape    = NULL,
                            scale    = NULL,
                            y        = NULL,
                            binf     = NULL,
                            bsup     = NULL,
                            initialize = function(type.prior=NA,opt.prior=NA,log=TRUE)
                            {
                              super$initialize(type.prior,opt.prior,log)
                              self$shape <- unlist(self$opt.prior)[1]
                              self$scale <- unlist(self$opt.prior)[2]
                              self$binf  <- 0
                              self$bsup  <- qgamma(0.9999,self$shape,1/self$scale)
                            },
                            prior = function(x=seq(0,self$bsup,length.out = 100))
                            {
                              if (self$log == FALSE)
                              {
                                return(1/(self$scale^self$shape*gamma(self$shape))*
                                         x^(self$shape-1)*exp(-(x/self$scale)))
                              } else
                              {
                                return(log(1/(self$scale^self$shape*gamma(self$shape))*
                                             x^(self$shape-1)*exp(-(x/self$scale))))
                              }
                            },
                            plot = function()
                            {
                              dplot <- data.frame(data=rgamma(n=100,shape=self$shape,scale=self$scale),type="prior")
                              p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                                geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                                theme_light()+xlab("")+ylab("")+
                                theme(legend.position=c(0.86,0.86),
                                      legend.text=element_text(face="bold",size = '20'),
                                      legend.title=element_blank(),
                                      legend.key=element_rect(colour=NA),
                                      axis.text=element_text(size=20))+
                                geom_hline(aes(yintercept = 0))
                              return(p)
                            }
                          ))



#Inverse gamma to be completed....
invGamma.class <- R6Class(classname = "invGamma.class",
                           inherit = prior.class,
                           public = list(
                             shape    = NULL,
                             scale    = NULL,
                             y        = NULL,
                             binf     = NULL,
                             bsup     = NULL,
                             initialize = function(type.prior=NA,opt.prior=NA,log=TRUE)
                             {
                               super$initialize(type.prior,opt.prior,log)
                               self$shape <- unlist(opt.prior)[1]
                               self$scale <- unlist(opt.prior)[2]
                               self$binf  <- 1/(self$shape*self$scale*(1-self$scale))
                               self$bsup  <- 1/(self$shape*self$scale*(1+self$scale))
                             },
                             prior = function(x=seq(0,1,length.out = 100))
                             {
                               if (self$log == FALSE)
                               {
                                 return(1/(1/(self$scale^self$shape*gamma(self$shape))*
                                          x^(self$shape-1)*exp(-(x/self$scale))))
                               } else
                               {
                                 return(log(1/(1/(self$scale^self$shape*gamma(self$shape))*
                                              x^(self$shape-1)*exp(-(x/self$scale)))))
                               }
                             },
                             plot = function()
                             {
                               dplot <- data.frame(data,1/ragamma(n=1000,self$shape,self$scale),type="prior")
                               p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                                 geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                                 theme_light()+xlab("")+ylab("")+ xlim(self$binf,self$bsup)+
                                 theme(legend.position=c(0.86,0.86),
                                       legend.text=element_text(face="bold",size = '20'),
                                       legend.title=element_blank(),
                                       legend.key=element_rect(colour=NA),
                                       axis.text=element_text(size=20))+
                                 geom_hline(aes(yintercept = 0))
                               return(p)
                             }
                           ))





