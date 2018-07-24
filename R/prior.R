#' A Reference Class to generate differents \code{\link{prior.class}} objects
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
#' @field opt.prior the characteristics of the selected prior
#' @export
prior.class <- R6Class(classname = "prior.class",
                           public = list(
                             type.prior  = NULL,
                             opt.prior   = NULL,
                             initialize = function(type.prior=NA,opt.prior=NA)
                             {
                               self$type.prior  <- type.prior
                               self$opt.prior   <- opt.prior
                               private$checkPrior()
                             },
                             centerText = function(x) {
                               width <- getOption("width")
                               out <- "=======================================\n"
                               l <- nchar(out)-2
                               l2 <- nchar(x)
                               l3 <- (l-l2)/2
                               temp <- rep("=",floor(l3))
                               mid <- paste(paste(temp[1:length(temp)],sep="",collapse = "")," ",x," ",
                                            paste(temp[1:length(temp)],sep="",collapse = ""),"\n",sep="")
                               ws <- rep(" ", floor((width - nchar(out))/2))
                               cat(ws, mid, sep = "")
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
                                initialize = function(type.prior=NA,opt.prior=NA)
                                {
                                  super$initialize(type.prior,opt.prior)
                                  self$mean  <- unlist(self$opt.prior)[1]
                                  self$var   <- unlist(self$opt.prior)[2]
                                  self$sd    <- sqrt(self$var)
                                  self$binf  <- self$mean - 5*sqrt(self$var)
                                  self$bsup  <- self$mean + 5*sqrt(self$var)
                                },
                                prior = function(x=seq(self$binf,self$bsup,length.out = 100))
                                {
                                  return(dnorm(x,mean=self$mean,sd=self$sd,log = TRUE))
                                },
                                plot = function()
                                {
                                  dplot <- data.frame(data=rnorm(n=1000,self$mean,sqrt(self$var)),type="prior")
                                  p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                                    geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                                    theme_light()+xlab("")+ylab("")+ xlim(c(self$binf,self$bsup))+
                                    theme(legend.title = element_blank(),legend.position = c(0.8,0.8),
                                          legend.background = element_rect(linetype="solid", colour ="grey"),
                                          axis.text=element_text(size=20))+
                                    geom_hline(aes(yintercept = 0))
                                  return(p)
                                },
                                print = function()
                                {
                                  cat("Call:\n")
                                  self$centerText("Gaussian prior")
                                  cat("\n")
                                  cat("Specifications: \n")
                                  print(data.frame(mean=self$mean,variance=self$var),row.names = FALSE)
                                }
                              ))


unif.class <- R6Class(classname = "unif.class",
                              inherit = prior.class,
                              public = list(
                                binf    = NULL,
                                bsup    = NULL,
                                y       = NULL,
                                initialize = function(type.prior=NA,opt.prior=NA)
                                {
                                  super$initialize(type.prior,opt.prior)
                                  self$binf <- unlist(self$opt.prior)[1]
                                  self$bsup <- unlist(self$opt.prior)[2]
                                },
                                prior = function(x=seq(self$binf,self$bsup,length.out = 100))
                                {
                                  return(dunif(x,min=self$binf,max=self$bsup,log=TRUE))
                                },
                                plot = function()
                                {
                                  xvals <- data.frame(data=c(self$binf,self$bsup),type="prior")
                                  ggplot(data.frame(xvals), aes(x = data,color=type)) +
                                    stat_function(fun = dunif, args = list(min =self$binf,max =self$bsup)
                                                  , geom = "area",fill = "red", alpha = 0.1) +
                                    theme_light()+xlab("")+ylab("")+ xlim(c(self$binf,self$bsup))+
                                    theme(legend.title = element_blank(),legend.position = c(0.8,0.8),
                                          legend.background = element_rect(linetype="solid", colour ="grey"),
                                          axis.text=element_text(size=20))+
                                    geom_hline(aes(yintercept = 0))
                                },
                                print = function()
                                {
                                  cat("Call:\n")
                                  self$centerText("Uniform prior")
                                  cat("\n")
                                  cat("Specifications: \n")
                                  print(data.frame(binf=self$binf,
                                                   bsup=self$bsup),row.names = FALSE)
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
                            initialize = function(type.prior=NA,opt.prior=NA)
                            {
                              super$initialize(type.prior,opt.prior)
                              self$shape <- unlist(self$opt.prior)[1]
                              self$scale <- unlist(self$opt.prior)[2]
                              self$binf  <- 0
                              self$bsup  <- qgamma(0.9999,self$shape,1/self$scale)
                            },
                            prior = function(x=seq(0,self$bsup,length.out = 100))
                            {
                              return(dgamma(x,scale=self$scale,shape=self$shape,log=TRUE))
                            },
                            plot = function()
                            {
                              dplot <- data.frame(data=rgamma(n=100,shape=self$shape,scale=self$scale),type="prior")
                              p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                                geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                                theme_light()+xlab("")+ylab("")+
                                xlim(c(self$binf,self$bsup))+
                                theme(legend.title = element_blank(),legend.position = c(0.8,0.8),legend.background = element_rect(
                                  linetype="solid", colour ="grey"),axis.text=element_text(size=20))+
                                geom_hline(aes(yintercept = 0))
                              return(p)
                            },
                            print = function()
                            {
                              cat("Call:\n")
                              self$centerText("Gamma prior")
                              cat("\n")
                              cat("Specifications: \n")
                              print(data.frame(scale=self$scale,shape=self$shape),row.names = FALSE)
                            }
                          ))





