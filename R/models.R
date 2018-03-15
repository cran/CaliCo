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
#' @field opt.emul a list of parameter for the surrogate \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (default value 1)}
#' \item{\strong{n.emul}}{ the number of points for constituting the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{binf}}{ the lower bound of the parameter vector (default value 0)}
#' \item{\strong{bsup}}{ the upper bound of the parameter vector (default value 1)}
#' \item{\strong{DOE}}{ design of experiments for the surrogate (default value NULL)}}
#' @field model the model choice (see \code{\link{model}} for more specification).
#' @field opt.disc a list of parameter for the discrepancy \itemize{
#' \item{\strong{kernel.type}}{ the kernel choosen for the Gaussian process}}
#'
#' @export
model.class <- R6Class(classname = "model.class",
                 public = list(
                   code      = NULL,
                   X         = NULL,
                   Yexp      = NULL,
                   n         = NULL,
                   d         = NULL,
                   binf      = NULL,
                   bsup      = NULL,
                   opt.emul = NULL,
                   model     = NULL,
                   opt.disc  = NULL,
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA,
                                         opt.emul=list(p=NA,n.emul=NA,type=NA,binf=NA,bsup=NA,DOE=NA))
                   {
                     self$code  <- code
                     self$X     <- X
                     self$Yexp  <- Yexp
                     self$n     <- length(Yexp)
                     if (is.null(dim(X)))
                     {
                       self$d   <- 1
                     } else
                     {
                       self$d   <- dim(X)[2]
                     }
                     self$opt.emul <- opt.emul
                     self$model     <- model
                     private$checkModels()
                     private$checkEmul()
                     private$checkCode()
                   }
                 ))

model.class$set("private","checkModels",
        function()
        {
          if (self$model != "model1" & self$model != "model2" & self$model != "model3" & self$model != "model4")
          {
            stop('Please elect a correct model')
          }
        })


model.class$set("private","checkEmul",
                function()
                  {
                  N <- c("p","n.emul","type","binf","bsup","DOE")
                  N2 <- names(self$opt.emul)
                  for (i in 1:length(N))
                  {
                    if(names(self$opt.emul)[i] != N[i])
                    {
                      stop(paste(N[i],"value is missing, please enter a correct value",sep=" "))
                    }
                  }
                })

model.class$set("private","checkCode",
                function()
                {
                  if (is.null(self$code))
                  {stop("Please enter a valid code")}
                })


model1.class <- R6Class(classname = "model1.class",
                        inherit = model.class,
                        public=list(
                        m.exp = NULL,
                        V.exp = NULL,
                        initialize=function(code=NA, X=NA, Yexp=NA, model=NA)
                        {
                          super$initialize(code, X, Yexp, model)
                        },
                        fun = function(theta,var)
                        {
                          y  <- self$code(self$X,theta)+rnorm(self$n,0,sqrt(var))
                          yc <- self$code(self$X,theta)
                          if (is.na(mean(yc)))
                          {stop('Wrong number of parameter in the function')}
                          return(list(y=y, yc=yc))
                        },
                        pred = function(theta,var,x.new)
                        {
                          if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                          y  <- self$code(x.new,theta)+rnorm(l,0,sqrt(var))
                          yc <- self$code(x.new,theta)
                          if (is.na(mean(yc)))
                          {stop('Wrong number of parameter in the function')}
                          return(list(y=y, yc=yc))
                        },
                        likelihood = function(theta,var)
                        {
                          self$m.exp = self$code(self$X,theta)
                          if (is.na(mean(self$m.exp)))
                          {stop('Wrong number of parameter in the function')}
                          self$V.exp = var*diag(self$n)
                          return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                        }
                        )
                        )


model1.class$set("public","plot",
                 function(theta,var,select.X=NULL)
                 {
                   if(is.null(dim(self$X))==FALSE & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   res <- self$fun(theta,var)
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){Xplot <- select.X}
                   binf = min(Xplot)
                   bsup = max(Xplot)
                   if (is.na(mean(res$yc)))
                   {stop('You have given the wrong number of parameter')}
                   gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),type="code result")
                   gg.data.noisy <- data.frame(y=res$y,x=seq(binf,bsup,length.out=length(res$yc)),
                                             type="noisy")
                   gg.data.exp  <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                              type="experiments")
                   gg.data <- rbind(gg.data,gg.data.noisy,gg.data.exp)
                   p <- ggplot(gg.data)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   return(p)
                 })

model1.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("No surrogate is selected")
                   cat("\n")
                   cat("No discrepancy is added")
                 }
)

model3.class <- R6Class(classname = "model3.class",
                        inherit = model1.class,
                        public=list(
                          funTemp  = NULL,
                          predTemp = NULL,
                          temp     = NULL,
                          disc     = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,
                                              opt.disc=list(kernel.type=NULL))
                          {
                            self$opt.disc  <- list(kernel.type=opt.disc$kernel.type)
                            if (is.null(self$opt.disc$kernel.type)==TRUE)
                            {
                              self$opt.disc$kernel.type="gauss"
                            }
                            super$initialize(code, X, Yexp, model)
                            self$funTemp  <- super$fun
                            self$predTemp <- super$pred
                          },
                          discrepancy = function(theta,thetaD,var,X=self$X)
                          {
                            y   <- self$funTemp(theta,var)$y
                            z   <- self$Yexp - y
                            Cov <- kernel.fun(X,thetaD[1],thetaD[2],self$opt.disc$kernel.type)
                            if (is.null(dim(X)) && length(X)==1)
                            {} else
                            {
                              p <- eigen(Cov)$vectors
                              e <- eigen(Cov)$values
                              if (all(e>0)){} else
                              {
                                e[which(e<0)] <- 1e-4
                              }
                              d <- diag(e)
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                Cov <- as.numeric(p)^2*d
                              } else
                              {
                                Cov <- t(p)%*%d%*%p
                              }
                            }
                            if (is.null(dim(X))){long <- length(X)}else{
                              long <- dim(X)[1]}
                            if (long==1)
                            {
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                biais <- rnorm(1,0,sqrt(Cov))
                              } else
                              {
                                biais <- rnorm(n=self$n,0,sqrt(Cov))
                              }
                            } else
                            {
                              biais <- mvrnorm(n=self$n,rep(0,long),Cov)
                              biais <- apply(biais,1,mean)
                            }
                            return(list(biais=biais,cov=Cov))
                          },
                          fun = function(theta,thetaD,var,X=self$X)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,X)
                            foo <- self$funTemp(theta,var)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
                          },
                          pred = function(theta,thetaD,var,x.new)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                            foo <- self$predTemp(theta,var,x.new)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
                          }
                          )
)


model3.class$set("public","plot",
                 function(theta,thetaD,var,select.X=NULL)
                 {
                   if(is.null(dim(self$X))==FALSE & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   res <- self$fun(theta,thetaD,var)
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){Xplot <- select.X}
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),type="code result")
                   gg.data.noisy <- data.frame(y=res$y,x=seq(binf,bsup,length.out=length(res$yc)),
                                               type="with discrepancy and noise")
                   gg.data.exp  <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                              type="experiments")
                   gg.data <- rbind(gg.data,gg.data.noisy,gg.data.exp)
                   p <- ggplot(gg.data)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   return(p)
                 })

model3.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("No surrogate is selected")
                   cat("\n\n")
                   cat("A discrepancy is added:\n")
                   cat(paste("Mean of the biais:",round(mean(self$disc$biais),5),"\n",sep=" "))
                   cat(paste("Covariance of the biais:",round(mean(self$disc$cov),3),"\n",sep=" "))
                   cat(paste("Kernel chossen: ",self$opt.disc$kernel.type,sep=""))
                 }
)



model3.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   self$m.exp <- self$code(self$X,theta)
                   temp <- self$fun(theta,thetaD,var)
                   self$V.exp <- var*diag(self$n) + temp$cov
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                              invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })



model2.class <- R6Class(classname = "model2.class",
                        inherit = model.class,
                        public = list(
                          n.emul   = NULL,
                          opt.emul = NULL,
                          GP       = NULL,
                          X        = NULL,
                          DOE      = NULL,
                          design   = NULL,
                          p        = NULL,
                          type     = NULL,
                          m.exp    = NULL,
                          V.exp    = NULL,
                        initialize = function(code=NA, X=NA, Yexp=NA, model=NA,opt.emul=NA)
                        {
                          super$initialize(code, X, Yexp, model,opt.emul)
                          self$opt.emul <- opt.emul
                          self$X        <- X
                          self$binf     <- opt.emul$binf
                          self$bsup     <- opt.emul$bsup
                          self$n.emul   <- opt.emul$n.emul
                          self$p        <- opt.emul$p
                          self$type     <- opt.emul$type
                          self$DOE      <- opt.emul$DOE
                          foo           <- self$surrogate()
                          self$GP       <- foo$GP
                          self$design   <- foo$doe
                          print("The surrogate has been set up, you can now use the function")
                        },
                        surrogate = function()
                        {
                          Xcr <- scale(self$X)
                          V   <- attr(Xcr,"scaled:scale")
                          M   <- attr(Xcr,"scaled:center")
                          Dim <- self$p+self$d
                          if (is.null(self$DOE)==FALSE)
                          {
                            D <- self$DOE
                            ### Generating the response
                            z <- matrix(nr=self$n.emul,nc=1)
                            for (i in 1:self$n.emul)
                            {
                              covariates <- as.matrix(D[i,1:(Dim-self$p)])
                              dim(covariates) <- c(1,self$d)
                              z[i] <- self$code(covariates,D[i,(Dim-self$p+1):Dim])
                            }
                            #z <- self$code(D[,1:(Dim-self$p)],D[,(Dim-self$p+1):Dim])
                            ### Converting D as a data.frame for the km function
                            D <- as.data.frame(D)
                            ### Creation of the Gaussian Process with estimation of hyperpameters
                            GP <- km(formula =~1, design=D, response = z,covtype = self$type)
                            return(list(GP=GP,doe=D))
                          }else
                          {
                            doe <- lhsDesign(self$n.emul,Dim)$design
                            doe <- maximinSA_LHS(doe)
                            doe <- doe$design
                            ### Getting back the value of the parameter generated by the DOE
                            binf.X <- apply(Xcr,2,min)
                            bsup.X <- apply(Xcr,2,max)
                            DOEtemp <- unscale(doe[,1:self$d],binf.X,bsup.X)
                            if (self$d==1)
                            {
                              for (i in 1:self$n.emul)
                              {
                                DOEtemp[i] <- DOEtemp[i]*V+M
                              }
                            } else {
                              for (i in 1:self$n.emul)
                              {
                                DOEtemp[i,] <- DOEtemp[i,]*V+M
                              }
                            }
                            if (length(self$binf)!=self$p | length(self$bsup)!=self$p)
                            {stop('Mismatch between the size of the upper, lower bounds and the parameter vector')}
                            doeParam <- unscale(doe[,(Dim-self$p+1):Dim],self$binf,self$bsup)
                            ### Matrix D contains the final value for the DOE
                            D <- cbind(DOEtemp,doeParam)
                            self$DOE <- D
                            ### Generating the response
                            z <- matrix(nr=self$n.emul,nc=1)
                            for (i in 1:self$n.emul)
                            {
                              covariates <- as.matrix(D[i,1:(Dim-self$p)])
                              dim(covariates) <- c(1,self$d)
                              z[i] <- self$code(covariates,D[i,(Dim-self$p+1):Dim])
                            }
                            #z <- self$code(D[,1:(Dim-self$p)],D[,(Dim-self$p+1):Dim])
                            ### Converting D as a data.frame for the km function
                            D <- as.data.frame(D)
                            #colnames(D) <- c("V1","V2","V3","V4","V5","V6")
                            ### Creation of the Gaussian Process with estimation of hyperpameters
                            GP <- km(formula =~1, design=D, response = z,covtype = self$type)
                            return(list(GP=GP,doe=D))
                          }
                        },
                        fun = function(theta,var)
                        {
                          if(self$p==1)
                          {
                            Xnew <- cbind(self$X,rep(theta,self$n))
                          } else
                          {
                            Xtemp <- matrix(rep(theta,rep(self$n,self$p)),nr=self$n,nc=self$p)
                            Xnew  <- cbind(self$X,Xtemp)
                          }
                          Xnew <- as.data.frame(Xnew)
                          names(Xnew) <- c("DOE","doeParam")
                          pr <- predict(self$GP,newdata=as.data.frame(Xnew),type="UK",
                                        cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          err <- rnorm(n=self$n,mean = 0,sd=sqrt(var))
                          return(list(y=pr$mean+err,Cov.GP=pr$cov,yc=pr$mean,lower=pr$lower95,upper=pr$upper95))
                        },
                        pred = function(theta,var,x.new)
                        {
                          if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                          if(self$p==1)
                          {
                            Xnew <- cbind(x.new,rep(theta,l))
                          } else
                          {
                            Xtemp <- matrix(rep(theta,rep(l,self$p)),nr=l,nc=self$p)
                            Xnew  <- cbind(x.new,Xtemp)
                          }
                          Xnew <- as.data.frame(Xnew)
                          #names(Xnew) <- c("DOE","doeParam")
                          pr <- predict(self$GP,newdata=as.data.frame(Xnew),type="UK",
                                        cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          if (length(pr$mean) == 1)
                          {
                            err <- rnorm(n=1,mean = 0,sd=sqrt(var))
                          } else
                          {
                            err <- rnorm(n=self$n,mean = 0,sd=sqrt(var))
                          }
                          return(list(y=pr$mean+err,Cov.GP=pr$cov,yc=pr$mean,lower=pr$lower95,upper=pr$upper95))
                        })
                        )


model2.class$set("public","likelihood",
                 function(theta,var)
                 {
                   temp <- self$fun(theta,var)
                   self$m.exp <- temp$yc
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   self$V.exp <- var*diag(self$n) + temp$Cov.GP
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                                invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })


model2.class$set("public","plot",
                 function(theta,var,points=FALSE,select.X=NULL)
                 {
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   if(self$d>1 & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   res <- self$fun(theta,var)
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){
                     Xplot <- select.X
                     if(points==TRUE)
                     {
                       points <- FALSE
                       print('The option point is automatically disabled because the representation would not have any sense')
                     }
                   }
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),
                                          lower=res$lower,upper=res$upper,type="Gaussian process",
                                         fill="90% credibility interval for the Gaussian process")
                   gg.data.exp <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),lower=res$lower,
                                            upper=res$upper,type="experiment",
                                            fill="90% credibility interval for the Gaussian process")
                   gg.data <- rbind(gg.data,gg.data.exp)
                   gg.points <- data.frame(x=self$DOE,y=self$code(self$DOE,theta))
                   p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     scale_fill_manual("",values=c("grey12"))+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   if (points==FALSE)
                   {
                     return(p)
                   } else
                   {
                     return(p+geom_jitter(data=gg.points,aes(x=x,y=y)))
                   }
                 })


model2.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("A surrogate had been set up:")
                   print(self$GP)
                   cat("\n")
                   cat("No discrepancy is added\n")
                 }
)



model4.class <- R6Class(classname = "model4.class",
                        inherit = model2.class,
                        public=list(
                          funC = NULL,
                          predTemp = NULL,
                          disc = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.emul=NA,
                                              opt.disc=list(kernel.type=NULL))
                          {
                            super$initialize(code, X, Yexp, model, opt.emul)
                            self$opt.disc  <- list(kernel.type=opt.disc$kernel.type)
                            if (is.null(self$opt.disc$kernel.type)==TRUE)
                            {
                              self$opt.disc$kernel.type="gauss"
                            }
                            self$funC     <- super$fun
                            self$predTemp <- super$pred
                          },
                          discrepancy = function(theta,thetaD,var,X=self$X)
                          {
                            y   <- self$funC(theta,var)$y
                            z   <- self$Yexp - y
                            Cov <- kernel.fun(X,thetaD[1],thetaD[2],self$opt.disc$kernel.type)
                            if (is.null(dim(X)) && length(X)==1)
                            {} else
                            {
                              p <- eigen(Cov)$vectors
                              e <- eigen(Cov)$values
                              if (all(e>0)){} else
                              {
                                e[which(e<0)] <- 1e-4
                              }
                              d <- diag(e)
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                Cov <- as.numeric(p)^2*d
                              } else
                              {
                                Cov <- t(p)%*%d%*%p
                              }
                            }
                            if (is.null(dim(X))){long <- length(X)}else{
                              long <- dim(X)[1]}
                            if (long==1)
                            {
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                biais <- rnorm(1,0,sqrt(Cov))
                              } else
                              {
                                biais <- rnorm(n=self$n,0,sqrt(Cov))
                              }
                            } else
                            {
                              biais <- mvrnorm(n=self$n,rep(0,long),Cov)
                              biais <- apply(biais,1,mean)
                            }
                            return(list(biais=biais,cov=Cov))
                          },
                          pred = function(theta,thetaD,var,x.new)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                            foo <- self$predTemp(theta,var,x.new)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
                          },
                          fun = function(theta,thetaD,var,X=self$X)
                          {
                            foo <- self$funC(theta,var)
                            self$disc <- self$discrepancy(theta,thetaD,var,X)
                            y <- foo$y
                            Cov.GP <- foo$Cov.GP
                            yc <- foo$yc
                            lower <- foo$lower
                            upper <- foo$upper
                            return(list(y=self$disc$biais+y,Cov.D=self$disc$cov,yc=yc,
                                        lower=lower,upper=upper,Cov.GP=Cov.GP))
                          })
)


model4.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   temp <- self$fun(theta,thetaD,var)
                   self$m.exp <- temp$yc
                   self$V.exp <- var*diag(self$n) + temp$Cov.GP +temp$Cov.D
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                                invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })



model4.class$set("public","plot",
                 function(theta,thetaD,var,points=FALSE,select.X=NULL)
                 {
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   if(self$d>1 & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   if(is.null(select.X))
                   {
                     res <- self$fun(theta,thetaD,var,self$X)
                   } else{
                     res <- self$fun(theta,thetaD,var,select.X)
                     }
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){Xplot <- select.X}
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),
                                         lower=res$lower,upper=res$upper,type="Gaussian process",
                                         fill="90% credibility interval for the Gaussian process")
                   gg.data.dis <- data.frame(y=res$y,x=seq(binf,bsup,length.out=length(res$yc)),lower=res$lower,
                                             upper=res$upper,type="Gaussian Process with discrepancy",
                                             fill="90% credibility interval for the Gaussian process")
                   gg.data.exp <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),lower=res$lower,
                                             upper=res$upper,type="Experiment",
                                             fill="90% credibility interval for the Gaussian process")
                   gg.data <- rbind(gg.data,gg.data.dis,gg.data.exp)
                   gg.points <- data.frame(x=self$DOE[,1],y=self$code(self$DOE,theta))
                   p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     scale_fill_manual("",values=c("grey12"))+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   if (points==FALSE)
                   {
                     return(p)
                   } else
                   {
                     return(p+geom_jitter(data=gg.points,aes(x=x,y=y)))
                   }
                 })

model4.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("A surrogate had been set up:")
                   print(self$GP)
                   cat("\n")
                   cat("A discrepancy is added:\n")
                   # cat(paste("Mean of the biais:",round(mean(self$disc$biais),5),"\n",sep=" "))
                   # cat(paste("Covariance of the biais:",round(mean(self$disc$cov),5),"\n",sep=" "))
                   cat("Chosen kernel:", self$opt.disc$kernel.type)
                 }
)

