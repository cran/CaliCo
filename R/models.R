#' A Reference Class to generates differents model objects
#'
#'
#' @description See the function \code{\link{model}} which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link{model}}... Other methods
#' should not be called as they are designed to be used during the calibration process.
#'
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field code a function which takes in entry X and theta
#' @field X the matrix of the forced variables
#' @field Yexp the experimental output
#' @field n the number of experiments
#' @field d the number of forced variables
#' @field binf the lower bound of the parameters for the DOE
#' @field bsup the upper bound of the parameters for the DOE
#' @field opt.gp a list of parameter for the surrogate (default NULL) \itemize{
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{DOE}}{ design of experiments for the surrogate (default value NULL)}}
#' @field opt.emul a list of parameter to establish the DOE (default NULL) \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (default value 1)}
#' \item{\strong{n.emul}}{ the number of points for constituting the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{binf}}{ the lower bound of the parameter vector (default value 0)}
#' \item{\strong{bsup}}{ the upper bound of the parameter vector (default value 1)}}
#' @field opt.sim a list of parameter containing output of the code and corresponding DOE \itemize{
#' \item{\strong{Ysim}}{ Output of the code}
#' \item{\strong{DOEsim}}{ DOE corresponding to the output of the code}}
#' @field model the model choice (see \code{\link{model}} for more specification).
#' @field opt.disc a list of parameter for the discrepancy \itemize{
#' \item{\strong{kernel.type}}{ the kernel chosen for the Gaussian process}}
#'
#' @export
model.class <- R6Class(classname = "model.class",
                 public = list(
                   code        = NULL, ## code variable
                   X           = NULL, ## forced variables
                   Yexp        = NULL, ## experiments
                   n           = NULL, ## length of the experiements
                   d           = NULL, ## number of forced variables
                   model       = NULL, ## statistical model elected
                   theta       = NULL, ## parameter vector
                   var         = NULL, ## variance of the measurement error
                   thetaD      = NULL, ## discrepancy parameter vector
                   n.cores     = NULL, ## number of computer cores
                   m.exp       = NULL, ## Mean for the likelihood
                   V.exp       = NULL, ## Variance for the likelihood
                   f.variables = NULL, ## variable indicating the presence
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA)
                   {
                     self$code    <- code
                     if (is.vector(X)) self$X <- matrix(X,ncol=1)
                     else self$X  <- as.matrix(X)
                     self$Yexp    <- Yexp
                     self$n       <- length(Yexp)
                     self$model   <- model
                     self$n.cores <- detectCores()
                     if (is.matrix(X)) {self$d <- ncol(X)} else{self$d <-1}
                     private$checkModels()
                     private$checkCode()
                     private$checkOptions()
                     private$testOnForcingVar()
                   },
                   active = list(
                     theta  = function(value) {return(value)},
                     thetaD = function(value) {return(value)},
                     var    = function(value) {return(value)}
                   ),
                   gglegend =function()
                   {
                     return(theme(legend.title = element_blank(),
                                  legend.position = c(0.2,0.8),
                                  legend.background = element_rect(linetype="solid", colour ="grey"),
                                  axis.title=element_text(size=20,family = "Helvetica",face = "italic"),
                                  axis.text = element_text(size=20),
                                  legend.text = element_text(size=10)))
                   }
                 ))

## Function that check if the right label of model is given
model.class$set("private","checkModels",
        function()
          ### Check if the chosen model is in the possible selection
        {
          if (!self$model %in% c("model1","model2","model3","model4"))
          {
            stop('Please elect a correct model',call. = FALSE)
          }
        })


## Check the content of the options
model.class$set("private","checkOptions",
                function()
                  ### Check if there is no missing in the options
                  {
                  test <- function(N,options)
                  {
                    N2 <- names(options)
                    for (i in 1:length(N))
                    {
                      if(names(options)[i] != N[i])
                      {
                        stop(paste(N[i],"value is missing, please enter a correct value",sep=" "),call. = FALSE)
                      }
                    }
                  }
                  if (self$model %in% c("model2","model4"))
                  {
                    if (!is.null(self$opt.emul)) test(c("p","n.emul","binf","bsup"),self$opt.emul)
                    test(c("type","DOE"),self$opt.gp)
                    if (!is.null(self$opt.sim)) test(c("Ysim","DOEsim"),self$opt.sim)
                  } else
                  {
                    if (self$model %in% c("model3","model4")){test(c("kernel.type"),self$opt.disc)}
                  }
                })

## Check if the code is present for Model1 or Model3
model.class$set("private","checkCode",
                function()
                  ### Check if the code is valid
                {
                  if (is.null(self$code) & self$model %in% c("model1","model3"))
                  {
                    stop("The code cannot be NULL if you chose model1 or model3",call. = FALSE)
                  }
                })

## Test of presence of frocing variables
model.class$set("private","testOnForcingVar",
                function()
                {
                  if (ncol(self$X) == 1 & nrow(self$X) == 1)
                  {
                    if (self$X == 0) self$f.variables <- 0 else self$f.variables <- 1
                  } else self$f.variables <- 1
                })

## Plot function for the models
model.class$set("public","plot",
                function(x,CI="all",...)
                {
                  if (missing(x)) stop("No x-axis selected, no graph is displayed",call. = FALSE)
                  if (is.matrix(x)){stop("please enter a correct x to plot your model",call. = FALSE)}
                  if (length(x)!= self$n){stop(paste("please enter a correct vector x of size",
                                                     self$n,sep=" "),call. = FALSE)}
                  df <- cbind(data.frame(y=self$Yexp,type="exp"),x=x)
                  if (self$model %in% c("model1","model2"))
                  {
                    if (!is.null(self$theta) & !is.null(self$var))
                    {
                      df2 <- cbind(self$model.fun(self$theta,self$var,X=self$X,CI),x=x)
                    } else
                      {
                        warning("no theta and var has been given to the model, experiments only are plotted",call.= FALSE)
                        df2 <- NULL
                      }
                  } else
                  {
                    if (!is.null(self$theta) & !is.null(self$thetaD) & !is.null(self$var))
                    {
                      df2 <- cbind(self$model.fun(self$theta,self$thetaD,self$var,X=self$X,CI),x=x)
                    } else
                      {
                        warning("no theta, thetaD and var has been given to the model, experiments only are plotted",call.= FALSE)
                        df2 <- NULL
                      }
                  }
                  if (!is.null(df2))
                  {
                    if (is.null(CI))
                    {
                      df <- rbind(df,df2)
                      p  <- ggplot(df)+geom_line(mapping = aes(x=x,y=y,color=type))+theme_light()+
                        xlab("")+ylab("")+self$gglegend() + scale_color_manual(values=c("red", "#000000"))
                    } else if (CI == "err")
                    {
                      if (self$model %in% c("model3","model4"))
                      {
                        df <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% discrepancy + noise")
                      } else
                      {
                        df <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% noise")
                      }
                      df2 <- df2[,names(df)]
                      df <- rbind(df,df2)
                      p  <- ggplot(df)+geom_ribbon(mapping = aes(x=x,ymin=q025,ymax=q975,fill=fill),alpha=0.4,linetype=1,
                                                   colour="skyblue3",size=0.5)+
                        geom_line(mapping = aes(x=x,y=y,color=type))+theme_light()+xlab("")+ylab("")+
                        scale_fill_manual(values = adjustcolor("skyblue3"))+
                        scale_color_manual(values=c("red", "#000000"))+self$gglegend()
                    } else if (CI == "GP")
                    {
                      if (self$model %in% c("model1","model3"))
                      {
                        stop("No Gaussian process used for the model1 and model2, no ggplot produced",call. = FALSE)
                      }
                      df <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% GP")
                      df <- rbind(df,df2)
                      p  <- ggplot(df)+geom_ribbon(mapping = aes(x=x,ymin=q025,ymax=q975,fill=fill),alpha=0.4,linetype=1,
                                                   colour="grey70",size=0.5)+
                        geom_line(mapping = aes(x=x,y=y,color=type))+theme_light()+xlab("")+ylab("")+
                        scale_fill_manual(values = adjustcolor("grey70"))+
                        scale_color_manual(values=c("red", "#000000"))+ self$gglegend()
                    } else if (CI == "all")
                    {
                      if (self$model %in% c("model1","model3"))
                      {
                        if (self$model == "model1"){df <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% noise")}
                        if (self$model == "model3")
                        {
                          df <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% discrepancy + noise")
                        }
                        df <- rbind(df,df2)
                        p  <- ggplot(df)+geom_ribbon(mapping = aes(x=x,ymin=q025,ymax=q975,fill=fill),alpha=0.4,linetype=1,
                                                     colour="skyblue3",size=0.5)+
                          geom_line(mapping = aes(x=x,y=y,color=type))+theme_light()+xlab("")+ylab("")+
                          scale_fill_manual(values = adjustcolor("skyblue3"))+
                          scale_color_manual(values=c("red", "#000000"))+self$gglegend()
                      } else
                      {
                        if (self$model == "model2")
                        {
                          df.gp  <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% GP")
                          df.n   <- cbind(df,q025=df2$q025n,q975=df2$q975n,fill="CI 95% noise")
                          df2.gp <- data.frame(y=df2$y,type=df2$type,x=x,q025=df2$q025,q975=df2$q975,fill="CI 95% GP")
                          df2.n  <- data.frame(y=df2$y,type=df2$type,x=x,q025=df2$q025n,q975=df2$q975n,
                                               fill="CI 95% noise")
                        } else
                        {
                          df.gp  <- cbind(df,q025=df2$q025,q975=df2$q975,fill="CI 95% GP")
                          df.n   <- cbind(df,q025=df2$q025n,q975=df2$q975n,fill="CI 95% discrepancy + noise")
                          df2.gp <- data.frame(y=df2$y,type=df2$type,x=x,q025=df2$q025,q975=df2$q975,fill="CI 95% GP")
                          df2.n  <- data.frame(y=df2$y,type=df2$type,x=x,q025=df2$q025n,q975=df2$q975n,
                                               fill="CI 95% discrepancy + noise")
                        }
                        df     <- rbind(df.n,df2.n)
                        df2    <- rbind(df.gp,df2.gp)
                        if (self$model == "model4")
                        {
                          col <- c("skyblue3","grey70")
                          Alpha <- c(0.8,0.3)
                        }else
                        {
                          col <- c("grey70","skyblue3")
                          Alpha <- c(0.3,0.8)
                        }
                        p <- ggplot(df)+
                          geom_ribbon(mapping = aes(x=x,ymin=q025,ymax=q975,fill=fill),
                                      alpha=0.8,linetype="twodash",colour="#999999",size=0.7)+
                          geom_ribbon(data=df2,mapping = aes(x=x,ymin=q025,ymax=q975,fill=fill),
                                      alpha=0.3,linetype="dotted",colour="#999999",size=0.7)+
                          geom_line(mapping = aes(x=x,y=y, color=type))+
                          scale_fill_manual(name = NULL,values = adjustcolor(col,alpha.f = 0.3))+
                          guides(fill = guide_legend(override.aes = list(alpha = Alpha)),
                                 colour= guide_legend(override.aes = list(colour = col)))+
                          scale_color_manual(values=c("red", "#000000"))+
                          xlab("")+ylab("")+theme_light()+self$gglegend()
                      }
                    } else
                    {
                      p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                        xlab("")+ylab("")+theme_light()+self$gglegend()+
                        scale_color_manual(values=c("red"))
                    }
                  } else
                  {
                    p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                      xlab("")+ylab("")+theme_light()+self$gglegend()+
                      scale_color_manual(values=c("red"))
                  }
                  return(p)
                })


## Print function
model.class$set("public","print",
                function(...)
                {
                  if (self$model %in% c("model3","model4"))
                  {
                    if (is.null(self$theta) | is.null(self$thetaD) | is.null(self$var)){
                      warning("\nThe discrepancy cannot be printed if no parameters are in the model",call. = FALSE)
                    } else
                    {
                      self$disc  <- self$discrepancy(self$theta,self$thetaD,self$var,self$X)
                      bias       <- summary(self$disc$bias)
                    }
                  }
                  cat("Call:\n")
                  print(self$model)
                  cat("\n")
                  cat("With the function:\n")
                  print(self$code)
                  cat("\n")
                  if (self$model %in% c("model1","model3"))
                  {
                    cat("No surrogate is selected")
                    cat("\n\n")
                    if (self$model == "model3")
                    {
                      cat("A discrepancy is added")
                      cat("\n\n")
                      cat("Summary of the bias mean:\n")
                      print(bias)
                      cat("Chosen kernel:", self$opt.disc$kernel.type)
                      cat(paste("\nCovariance of the bias:",round(mean(self$disc$cov),3),"\n\n",sep=" "))
                    } else {
                      cat("No discrepancy is added")
                    }
                  } else
                  {
                    cat("A surrogate had been set up:")
                    if (self$f.variables == 0)
                    {
                      if (self$lenCode > 1)
                      {
                        print("In time series modeling, a Gaussian process is available for each time step. You can access each one of them by yourmodel$GP")
                      } else
                      {
                        print(self$GP)
                      }
                    } else
                    {
                      print(self$GP)
                    }
                    cat("\n")
                    if (self$model == "model2")
                    {
                      cat("No discrepancy is added")
                    } else
                    {
                      if (is.null(self$theta) | is.null(self$thetaD) | is.null(self$var)){
                        cat("Summary of the bias mean:\n")
                        cat("No discrepancy in the print function")
                      } else
                      {
                        cat("Summary of the bias mean:\n")
                        print(bias)
                        cat("Chosen kernel:", self$opt.disc$kernel.type)
                        cat(paste("\nCovariance of the bias:",round(mean(self$disc$cov),3),"\n\n",sep=" "))
                      }
                    }
                  }
                })

## Discrepancy function
model.class$set("public","discrepancy",
                ## Define method that generates a discrepancy
                function(theta,thetaD,var,X=self$X)
                {
                  if (self$model=="model3")
                  {
                    y <- self$model1.fun(theta,var,X)$y
                  } else y <- self$model2.fun(theta,var,X)$y
                  z   <- self$Yexp - y
                  ## Compute the discrepancy covariance
                  Cov <- kernel.fun(X,thetaD[1],thetaD[2],self$opt.disc$kernel.type)
                  if (is.vector(X) & length(X)==1)
                  {} else
                  {
                    p <- eigen(Cov)$vectors
                    e <- eigen(Cov)$values
                    if (all(e>0)){} else
                    {
                      e[which(e<0)] <- .Machine$double.eps
                    }
                    d <- diag(e)
                    if (nrow(p) == 1 & ncol(p) == 1)
                    {
                      Cov <- as.numeric(p)^2*d
                    } else
                    {
                      Cov <- p%*%d%*%t(p)
                    }
                  }
                  if (is.vector(X)){long <- length(X)}else{
                    long <- nrow(X)}
                  if (long==1)
                  {
                    if (nrow(p) == 1 & ncol(p) == 1)
                    {
                      bias <- rnorm(1,0,sqrt(Cov))
                    } else
                    {
                      bias <- rnorm(n=self$n,0,sqrt(Cov))
                    }
                  } else
                  {
                    bias <- mvrnorm(10000,rep(0,long),Cov)
                    dim(bias) <- c(long,10000)
                    qq <- apply(bias,1,quantile,c(0.025,0.5,0.975))
                    bias <- qq[2,]
                    lower <- qq[1,]
                    upper <- qq[3,]
                  }
                  return(list(bias=bias,cov=Cov,lower=lower,upper=upper))
                })


################################## Model 1 definition #######################################
model1.class <- R6Class(classname = "model1.class",
                        inherit = model.class,
                        public=list(
                        initialize=function(code=NA, X=NA, Yexp=NA, model=NA)
                        {
                          ### Initialize from model.class
                          super$initialize(code, X, Yexp, model)
                        },
                        model.fun = function(theta,var,X=self$X,CI="err")
                        {
                          ### Function that generates the output of the model. If CI=TRUE, it computes the credibility
                          ### intervals of the white Gaussian noise
                          #y <- matrix(nr=100,nc=self$n)
                          # for(i in 1:100){y[i,] <- self$code(X,theta)+rnorm(self$n,0,sqrt(var))}
                          # qq <- apply(y,2,quantile,c(0.05,0.5,0.95))
                          y <- self$code(X,theta)
                          qq025 <- y - 2*sqrt(var)
                          qq975 <- y + 2*sqrt(var)
                          if (is.null(CI))
                          {
                            # df <- data.frame(y=qq[2,],type="model output")
                            df <- data.frame(y=y,type="model output")
                          } else if (CI=="err" | CI == "all")
                          {
                            # df <- data.frame(y=qq[2,],type="model output",q025=qq[1,],q975=qq[3,],fill="CI 95% noise")
                            df <- data.frame(y=y,type="model output",q025=qq025,q975=qq975,fill="CI 95% noise")
                          } else
                          {
                            warning("The argument for the credibility interval is not valid and no credibility interval will be displayed",call. = FALSE)
                            # df <- data.frame(y=qq[2,],type="model output")
                            df <- data.frame(y=y,type="model output")
                          }
                          return(df)
                        }
                        # prediction.fun = function(theta,var,x.new)
                        # {
                        #   ### Prediction function is the function to use when applying on a new data set
                        #   if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                        #   y  <- self$code(x.new,theta)+rnorm(l,0,sqrt(var))
                        #   df <- data.frame(pr=y,type="predicted")
                        #   return(df)
                        #}
                        )
)

## likelihood function
model1.class$set("public","likelihood",
                function(theta,var)
                {
                  ### Log-Likelihood
                  self$m.exp <- self$code(self$X,as.vector(theta))
                  self$V.exp <- var*diag(self$n)
                  return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                         -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                })



##################################### Model 3 definition ##################################

## model main functions
model3.class <- R6Class(classname = "model3.class",
                        inherit = model1.class,
                        public=list(
                          model1.fun        = NULL, ## model function from model1
                          model1.prediction = NULL, ## prediction function from model1
                          opt.disc          = NULL, ## discrepancy options
                          disc              = NULL, ## discrepancy field
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.disc=list(kernel.type=NULL))
                          {
                            ## Check if the opt.emul option is filled if it is not a gaussian kernel is picked
                            if (is.null(opt.disc$kernel.type))
                            {
                              warning("default value is selected. The discrepancy will have a gauss covariance structure",call. = FALSE)
                              self$opt.disc$kernel.type="gauss"
                            } else
                            {
                              self$opt.disc  <- opt.disc
                            }
                            ## Initialize with the model1.class which initialize from model.class
                            super$initialize(code, X, Yexp, model)
                            ## Store the model1 functions
                            self$model1.fun        <- super$model.fun
                          },
                          model.fun = function(theta,thetaD,var,X=self$X,CI="err")
                          {
                            res.model1 <- self$model1.fun(theta,var,X,CI=CI)
                            self$disc  <- self$discrepancy(theta,thetaD,var,X)
                            if (is.null(CI))
                            {
                              df <- data.frame(y=res.model1$y,type="model output")
                            } else if (CI=="err" | CI == "all")
                            {
                              # qq025 <- self$disc$lower + res.model1$q025
                              # qq975 <- self$disc$upper + res.model1$q975
                              qq025 <- res.model1$y - 2*sqrt(self$thetaD[1] + self$var)
                              qq975 <- res.model1$y + 2*sqrt(self$thetaD[1] + self$var)
                              df <- data.frame(y=res.model1$y,type="model output",q025=qq025,q975=qq975,
                                               fill="CI 95% discrepancy + noise")
                            } else
                            {
                              df <- 0
                            }
                            return(df)
                          }
                          # prediction.fun = function(theta,thetaD,var,x.new)
                          # {
                          #   self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                          #   foo <- self$model1.fun(theta,var,x.new)
                          #   y <- foo$y
                          #   df <- data.frame(y=self$disc$bias+y,cov=self$disc$cov,type="predicted")
                          #   return(df)
                          # }
                          )
)

## likelihood function
model3.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   # Log-Likelihood
                   self$disc  <- self$discrepancy(theta,thetaD,var,X=self$X)
                   self$m.exp <- self$code(self$X,as.vector(theta))
                   self$V.exp <- var*diag(self$n) + self$disc$cov
                   return(as.numeric(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                          -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp)))
                 })



################################## Model 2 definition ####################################

## model 2 main functions
model2.class <- R6Class(classname = "model2.class",
                        inherit = model.class,
                        public = list(
                          opt.emul = NULL, ## DOE creation options
                          opt.sim  = NULL, ## options if the user possess the design and the output
                          opt.gp   = NULL, ## GP options
                          case     = NULL, ## case wanted by the user (depending on options)
                          doe      = NULL, ## DOE used for the surrogate
                          z        = NULL, ## output of the code for the DOE
                          GP       = NULL, ## current Gaussian process emulated
                          p        = NULL, ## number of parameters
                          lenCode  = NULL, ## length of code output
                        initialize = function(code=NA, X=NA, Yexp=NA, model=NA,opt.gp=NULL,opt.emul=NULL,
                                              opt.sim=NULL,...)
                        {
                          if (missing(opt.sim) | is.null(opt.sim)) self$case <- 1
                          if ((missing(opt.emul) | is.null(opt.emul)) & (missing(opt.sim) | is.null(opt.sim)))
                            self$case <- 2
                          if ((missing(opt.emul) | is.null(opt.emul)) & !(missing(opt.sim) | is.null(opt.sim)))
                            self$case <- 3
                          self$opt.gp   <- opt.gp
                          self$opt.emul <- opt.emul
                          self$opt.sim  <- opt.sim
                          ## initialize from model.class
                          super$initialize(code, X, Yexp, model)
                          if (is.null(self$code))
                          {
                            if (is.null(opt.sim))
                            {
                              print("The numerical code is desabled, please fill the opt.sim option")
                            }
                          }
                          self$GP     <- self$surrogate()
                          print("The surrogate has been set up, you can now use the function")
                        },
                        ## model function
                        model.fun = function(theta,var,X=self$X,CI="all")
                        {
                          X  <- as.matrix(X)
                          self$p <- length(theta)
                          if (ncol(X) != self$d){stop("please enter a correct X",call. = FALSE)}
                          if (self$p == 1){
                            D <- cbind(X,rep(theta,nrow(X)))
                          } else
                          {
                            D  <- cbind(X,matrix(rep(theta,rep(nrow(X),self$p)),
                                              nr=nrow(X),nc=self$p))
                          }
                          if (self$f.variables == 0)
                          {
                            D <- D[,-1]
                            if (self$lenCode>1)
                            {
                              pr <- list()
                              for (i in 1:length(self$GP))
                              {
                                prTemp <- predict(self$GP[[i]],newdata=as.data.frame(t(D)),type="UK",
                                                   cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                                pr$mean[i] <- prTemp$mean
                                pr$lower95[i] <- prTemp$lower95
                                pr$upper95[i] <- prTemp$upper95
                                pr$cov[i] <- prTemp$cov
                              }
                            } else pr <- predict(self$GP,newdata=as.data.frame(D),type="UK",
                                          cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          } else pr <- predict(self$GP,newdata=as.data.frame(D),type="UK",
                                               cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          # nugget <- mvrnorm(n=100,pr$mean,diag(var,length(pr$mean)))
                          qq025 <- pr$mean -2*sqrt(var)
                          qq975 <- pr$mean +2*sqrt(var)
                          if (is.null(CI))
                          {
                            # qq <- apply(nugget,2,mean)
                            # df <- data.frame(y=qq,type="model output")
                            df <- data.frame(y=pr$mean,type="model output")
                          } else if (CI == "all" | CI == "err")
                          {
                            # qq <- apply(nugget,2,quantile,c(0.05,0.5,0.95))
                            if (CI == "all")
                            {
                              # df  <- data.frame(y=qq[2,],type="model output",q025n=qq[1,],q975n=qq[3,],
                              #                   q025=pr$lower95,q975=pr$upper95)
                              df  <- data.frame(y=pr$mean,type="model output",q025n=qq025,q975n=qq975,
                                                q025=pr$lower95,q975=pr$upper95)
                            } else
                            {
                              # df  <- data.frame(y=qq[2,],type="model output",q025=qq[1,],q975=qq[3,],
                              #                   fill="CI 95% noise")
                              df  <- data.frame(y=pr$mean,type="model output",q025=qq025,q975=qq975,
                                                fill="CI 95% noise")
                            }
                          } else if (CI == "GP")
                          {
                            # qq <- apply(nugget,2,quantile,c(0.05,0.5,0.95))
                            # df  <- data.frame(y=qq[2,],type="model output",q025=pr$lower95,
                            #                   q975=pr$upper95, fill="CI 95% GP")
                            df  <- data.frame(y=pr$mean,type="model output",q025=pr$lower95,
                                              q975=pr$upper95, fill="CI 95% GP")
                          } else
                          {
                            warning("The argument for the credibility interval is not valid and no credibility interval will be displayed",call. = FALSE)
                            df <- 0
                          }
                          return(df)
                        })
                        )

model2.class$set("public","surrogate",
                 function()
                 {
                   if (self$case == "1" | self$case =="2")
                   {
                     if (self$case == "1")
                     {
                       ## Dim is the dimension of H*Theta
                       if (self$f.variables == 0)
                       {
                         Dim <- self$opt.emul$p
                         ## Creation of the maximin LHS
                         doe <- lhsDesign(self$opt.emul$n.emul,Dim)$design
                         doe <- maximinSA_LHS(doe)$design
                         ## Going back to the original space of Theta
                         doe.theta <- unscale(doe,self$opt.emul$binf,
                                              self$opt.emul$bsup)
                         ## Generate the doe
                         self$doe  <- doe.theta
                       } else
                       {
                         Dim <- self$opt.emul$p+self$d
                         ## Creation of the maximin LHS
                         doe <- lhsDesign(self$opt.emul$n.emul,Dim)$design
                         doe <- maximinSA_LHS(doe)$design
                         ## Get the boundaries of X
                         binf.X <- apply(self$X,2,min)
                         bsup.X <- apply(self$X,2,max)
                         ## Going back to the original space H and Theta
                         doe.X <- unscale(doe[,c(1:self$d)],binf.X,bsup.X)
                         doe.theta <- unscale(doe[,c((self$d+1):Dim)],self$opt.emul$binf,
                                              self$opt.emul$bsup)
                         ## Generate the doe
                         self$doe  <- cbind(doe.X,doe.theta)
                         self$p <- self$opt.emul$p
                       }
                     } else
                     {
                       ###
                       self$doe <- self$opt.gp$DOE
                       Dim <- ncol(self$doe)
                       if (self$f.variables == 0){self$p <- Dim}
                       else {self$p <- Dim - self$d}
                     }
                     ## Compute the output of the code for the doe
                     self$z <- NULL
                     if (self$f.variables == 0)
                     {
                       if (is.matrix(self$doe))
                       {
                         self$lenCode <- length(self$code(0,self$doe[1,]))
                         if (self$lenCode>1) self$z <- matrix(nr=nrow(self$doe),nc=self$lenCode)
                         l <- nrow(self$doe)
                       } else
                       {
                         self$lenCode <- length(self$code(0,self$doe))
                         if (self$lenCode>1) self$z <- matrix(nr=length(self$doe),nc=self$lenCode)
                         l <- length(self$doe)
                       }
                     } else
                     {
                       l <- nrow(self$doe)
                     }
                     for (i in 1:l)
                     {
                       if (!self$f.variables == 0)
                       {
                         covariates <- as.matrix(self$doe[i,1:(Dim-self$p)])
                         if (self$d == 1) {} else dim(covariates) <- c(1,self$d)
                         self$z <- c(self$z,self$code(covariates,self$doe[i,(Dim-self$p+1):Dim]))
                       } else if (self$lenCode >1)
                       {
                         if (is.matrix(self$doe)) self$z[i,] <- self$code(0,self$doe[i,]) else
                           self$z[i,] <- self$code(0,self$doe[i])
                       } else
                       {
                         self$z <- c(self$z,self$code(0,self$doe[i,]))
                       }
                     }
                   } else if (self$case =="3")
                   {
                      if (self$f.variables == 0){self$p <- ncol(self$opt.sim$DOEsim)}
                      else {self$p <- ncol(self$opt.sim$DOEsim)-self$d}
                      self$doe <- self$opt.sim$DOEsim
                      self$z <- self$opt.sim$Ysim
                      if (is.matrix(self$z)) self$lenCode <- ncol(self$z)
                   }
                   ## Create the Gaussian Process corresponding
                   if (!self$f.variables == 0)
                   {
                      GP <- km(formula =~1, design=self$doe, response = self$z,covtype = self$opt.gp$type)
                   } else if (self$lenCode > 1)
                   {
                      GP <- list()
                      for (i in 1:self$lenCode){
                          GP[[i]] <- km(formula =~1, design=as.data.frame(self$doe), response = self$z[,i],
                                        covtype = self$opt.gp$type)
                   }
                   } else
                   {
                      GP <- km(formula =~1, design=self$doe, response = self$z,covtype = self$opt.gp$type)
                   }
                   return(GP)
                 })


## Likelihood
model2.class$set("public","likelihood",
                 function(theta,var)
                 {
                   D  <- cbind(self$X,matrix(rep(theta,rep(nrow(self$X),self$p)),
                                        nr=nrow(self$X),nc=self$p))
                   pr <- predict(self$GP,newdata=as.data.frame(D),type="UK",
                                 cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                   self$m.exp <- pr$mean
                   self$V.exp <- var*diag(self$n) + pr$cov
                   return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                          -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(Yexp-self$m.exp))
                 })



################################## Model 4 definition ####################################

## model 4 main functions
model4.class <- R6Class(classname = "model4.class",
                        inherit = model2.class,
                        public=list(
                          model2.fun = NULL, ## function from model2
                          opt.disc   = NULL, ## discrepancy options
                          disc       = NULL, ## discrepancy field
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.disc=NULL,opt.gp=NULL,
                                              opt.emul=NULL,opt.sim=NULL,...)
                          {
                            if (missing(opt.sim)) self$case <- 1
                            if (missing(opt.emul) & missing(opt.sim)) self$case <- 2
                            if (missing(opt.emul) & !missing(opt.sim)) self$case <- 3
                            self$opt.gp   <- opt.gp
                            self$opt.emul <- opt.emul
                            self$opt.sim  <- opt.sim
                            self$opt.disc <- opt.disc
                            super$initialize(code, X, Yexp, model,opt.gp=self$opt.gp, opt.emul=self$opt.emul,
                                             opt.sim=self$opt.sim)
                            ## Check if the opt.emul option is filled if it is not a gaussian kernel is picked
                            if (is.null(opt.disc$kernel.type)==TRUE)
                            {
                              warning("default value is selected. The discrepancy will have a gaussian covariance structure",call.=FALSE)
                              self$opt.disc$kernel.type="gauss"
                            } else
                            {
                              self$opt.disc <- opt.disc
                            }
                            self$model2.fun <- super$model.fun
                          },
                          model.fun = function(theta,thetaD,var,X=self$X,CI="all")
                          {
                            res.model2 <- self$model2.fun(theta,var)
                            self$disc <- self$discrepancy(theta,thetaD,var,X)
                            if (is.null(CI))
                            {
                              df <- data.frame(y=res.model2$y,type="model output")
                            } else if (CI == "err")
                            {
                              qq025 <- res.model2$y - 2*sqrt(self$thetaD[1] + self$var)
                              qq975 <- res.model2$y + 2*sqrt(self$thetaD[1] + self$var)
                              df    <- data.frame(y=res.model2$y,type="model output",q025=qq025,q975=qq975,
                                               fill="CI 95% discrepancy + noise")
                            } else if (CI == "GP")
                            {
                              df    <- data.frame(y=res.model2$y,type="model ouput",q025=res.model2$q025,
                                                      q975 = res.model2$q975,fill="CI 95% GP")
                            } else if (CI == "all")
                            {
                              qq025 <- res.model2$y - 2*sqrt(self$thetaD[1] + self$var)
                              qq975 <- res.model2$y + 2*sqrt(self$thetaD[1] + self$var)
                              df    <- data.frame(y=res.model2$y,type="model ouput",q025=res.model2$q025,
                                                      q975 = res.model2$q975,fill="CI 95% GP",q025n=qq025,
                                                      q975n=qq975)
                            } else
                            {
                              warning("The argument for the credibility interval is not valid and no credibility interval will be displayed",call. = FALSE)
                              df <- 0
                            }
                            return(df)
                          })
)

## Likelihood
model4.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   D  <- cbind(self$X,matrix(rep(theta,rep(nrow(self$X),self$p)),
                                             nr=nrow(self$X),nc=self$p))
                   pr <- predict(self$GP,newdata=as.data.frame(D),type="UK",
                                 cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                   dd <- self$discrepancy(theta,thetaD,var,self$X)
                   self$m.exp <- pr$mean
                   self$V.exp <- var*diag(self$n) + pr$cov +dd$cov
                   return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                          -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })
