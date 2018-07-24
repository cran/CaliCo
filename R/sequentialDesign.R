#' A Reference Class to generate a better Model2 or Model4 \code{\link{model.class}} objects
#'
#' @description This class generates a new \code{\link{model.class}} for Model4 and Model2. Based on the previous
#' estimation of the Gaussian process in the function \code{\link{model}}, the design of experiments previously used
#' is improved according to [Damblin et al. 2018]. The aim is to reduce the
#' error produced by the initial estimation of the Gaussian process by fortifying the initial DOE. The method consists
#' in proposing new points based on the expectancy improvement criterion.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field doe.init the initial DOE used to fit the first Gaussian process
#' @field GP.init the initial Gaussian process generated in \code{\link{model}} function
#' @field GP.new the new Gaussian process fortified with the new design points
#' @field p the number of parameter
#' @field md the initial model
#' @field md.new the new model
#' @field mdfit the initial calibrated model
#' @field mdfit.new the new calibrated model
#' @field X the data set
#' @field m minimum of the sum of squares used in the algorithm
#' @references DAMBLIN, Guillaume, BARBILLON, Pierre, KELLER, Merlin, et al. Adaptive numerical designs for the
#'  calibration of computer codes. SIAM/ASA Journal on Uncertainty Quantification, 2018, vol. 6, no 1, p. 151-179.
#' @export
seqDesign.class <- R6Class(classname = "seqDesign.class",
                       public = list(
                         doe.init  = NULL,
                         GP.init   = NULL,
                         doe.new   = NULL,
                         GP.new    = NULL,
                         p         = NULL,
                         md        = NULL,
                         md.new    = NULL,
                         mdfit     = NULL,
                         mdfit.new = NULL,
                         optimGrid = NULL,
                         X         = NULL,
                         m         = NULL,
                         initialize = function(md = NA,pr = NA,opt.estim = NA,k = NA)
                         {
                           if (!md$model %in% c("model2","model4"))
                           {
                             stop("The sequential design is available only for mode 2 and 4", call. = FALSE)
                           }
                           if (is.null(md$code))
                           {
                             stop("The sequential design is available only if the code field is given", call.=FALSE)
                           }
                           ### Uniquement pour model2
                           #### get the DOE and the GP from the model
                           self$doe.init  <- md$doe
                           self$GP.init   <- md$GP
                           self$md        <- md
                           self$p         <- ncol(self$doe.init) - self$md$d
                           X              <- md$X
                           binf           <- apply(self$doe.init,2,min)[(self$md$d+1):(self$md$d+self$p)]
                           bsup           <- apply(self$doe.init,2,max)[(self$md$d+1):(self$md$d+self$p)]
                           self$mdfit     <- calibrate(md,pr,opt.estim)
                           ### get the maximum a posteriori
                           thetaHat       <- estimators(self$mdfit)$MAP[1:self$p]
                           ### find the x* from Crit function
                           y.new          <- self$Crit(thetaHat)
                           self$doe.new   <- rbind(self$doe.init,y.new)
                           ## Update the new DOE GP with the new doe
                           z.new          <- c(md$z,md$code(as.matrix(t(y.new[1:self$md$d])),
                                                           as.matrix(t(y.new[(self$md$d+1):(self$md$d+self$p)]))))
                           self$GP.new    <- km(formula =~1, design=self$doe.new,
                                               response = z.new,covtype = md$opt.gp$type)
                           ## compute the sum of squares
                           self$m         <- self$SS(thetaHat)
                           ## enter the loop
                           thetaTemp      <- thetaHat
                           for (i in 1:k)
                           {
                             self$optimGrid <- unscale(maximinSA_LHS(lhsDesign(100*self$p,self$p)$design)$design,binf,bsup)
                             EIparallel <- function(i)
                             {
                               return(self$EI(self$optimGrid[i,]))
                             }
                             if (Sys.info()[['sysname']]=="Windows")
                             {
                               n.cores <- 1
                             } else
                             {
                               n.cores <- 2
                             }
                             EIvector  <- unlist(mclapply(c(1:(100*self$p)),EIparallel,mc.cores = n.cores))
                             if (all(EIvector == 0)) break
                             thetaHatNew  <- self$optimGrid[which(EIvector==max(EIvector)),]
                             if (all(thetaHatNew == thetaTemp))
                             {
                               break
                             } else
                             {
                               thetaTemp    <- thetaHatNew
                               y.new        <- self$Crit(thetaHatNew)
                               self$doe.new <- rbind(self$doe.new,y.new)
                               z.new        <- c(z.new,md$code(as.matrix(t(y.new[1:self$md$d])),
                                                               as.matrix(t(y.new[(self$md$d+1):(self$md$d+self$p)]))))
                               self$GP.new  <- km(formula =~1, design=self$doe.new,
                                                  response = z.new,covtype = md$opt.gp$type)
                               self$m       <- c(self$m[1:k],min(c(self$m[1:k],self$SS(thetaHatNew))))
                             }
                           }
                           cat("\n")
                           cat("======================================\n")
                           print(paste(i," new points generated",sep=""))
                           cat("\n======================================")
                           cat("\n")
                           if (self$md$model == "model2")
                           {
                             self$md.new    <- model(code=NULL,X = X,Yexp = self$md$Yexp,model = self$md$model,
                                                opt.gp=self$md$opt.gp, opt.sim=list(Ysim=z.new,DOEsim=self$doe.new))
                           } else
                           {
                             self$md.new    <- model(code=NULL,X = X,Yexp = self$md$Yexp,model = self$md$model,
                                                     opt.gp=self$md$opt.gp, opt.sim=list(Ysim=z.new,DOEsim=self$doe.new),
                                                     opt.disc = self$md$opt.disc)
                           }
                           self$mdfit.new <- calibrate(md = self$md.new,pr=pr, opt.estim = opt.estim)
                           invisible(list(md=self$md.new,mdfit=self$mdfit.new))
                         },
                         Crit = function(theta)
                         {
                           cov <- numeric(nrow(self$md$X))
                           for (i in 1:nrow(self$md$X))
                           {
                             newD <- c(self$md$X[i,],theta)
                             pr   <- predict(self$GP.init,newdata=as.data.frame(t(newD)),type="UK",
                                             interval="confidence",checkNames=FALSE)
                             cov[i]  <- pr$sd
                           }
                           return(c(self$md$X[which(cov==max(cov)),],theta))
                         },
                         fun.new = function(X,theta)
                         {
                           if (self$md$d == 1 & self$p == 1)
                           {
                             new.design <- cbind(X,rep(theta,length(X)))
                           } else if (self$md$d != 1 & self$p == 1)
                           {
                             new.design <- cbind(X,rep(theta,nrow(X)))
                           } else
                           {
                             new.design <- cbind(X,t(replicate(nrow(X),theta)))
                           }
                           pr           <- predict(self$GP.new, newdata=new.design,type="UK",
                                                   checkNames=FALSE,cov.compute=TRUE)
                           return(pr)
                         },
                         SS = function(theta)
                         {
                           y <- self$fun.new(self$md$X,theta)$mean
                           return(sum((self$md$Yexp-y)^2))
                         },
                         EI = function(theta, k=1)
                         {
                           pr   <- self$fun.new(self$md$X,theta)
                           S    <- mvrnorm(10000,pr$mean,pr$cov)
                           y    <- matrix(self$md$Yexp,nrow=10000,ncol = length(self$md$Yexp),byrow = TRUE)
                           SSsim <- rowSums((y-S)^2)
                           return(mean(apply(cbind(self$m[k]-SSsim,0),1,max)))
                         }
                       ))

seqDesign.class$set("public","plot",
                    function(x,graph="all",...)
                    {
                      res <- list()
                      if (graph == FALSE)
                      {
                        if (self$p != 1)
                        {
                          n        <- combinations(n=self$p,r=2,repeats=FALSE)
                          p        <- list()
                          for (i in 1:nrow(n))
                          {
                            label1 <- n[i,1]
                            label2 <- n[i,2]
                            ind1   <- n[i,1]+self$md$d
                            ind2   <- n[i,2]+self$md$d
                            df     <- data.frame(x=self$doe.init[,ind1],y=self$doe.init[,ind2],col="Initial DOE")
                            df2    <- data.frame(x=self$doe.new[-c(1:nrow(self$doe.init)),ind1],
                                                 y=self$doe.new[-c(1:nrow(self$doe.init)),ind2],col="Additional points")
                            df     <- rbind(df,df2)
                            p[[i]] <- ggplot(df, aes(x = x,y = y, color=col)) + geom_point(shape=3) + theme_light() +
                              xlab(substitute(theta[label1])) + ylab(substitute(theta[label2])) +
                              theme(legend.title = element_blank(),legend.position = c(0.2,0.8),
                                    legend.background = element_rect(linetype="solid", colour ="grey"))
                          }
                          res$doe <- p
                          t1     <- plot(self$mdfit,x,graph=NULL)$out + ggtitle("Before sequential design")
                          t2     <- plot(self$mdfit.new,x,graph=NULL)$out + ggtitle("After sequential design")
                          t3     <- plot(self$mdfit,x,graph=NULL)$dens
                          t4     <- plot(self$mdfit.new,x,graph=NULL)$dens
                          res$res <- list(dens=list(t3,t4),results=list(t1,t2))
                        } else
                        {
                          t1       <- plot(self$mdfit,x,graph=NULL)$out + ggtitle("Before sequential design")
                          t2       <- plot(self$mdfit.new,x,graph=NULL)$out + ggtitle("After sequential design")
                          t3     <- plot(self$mdfit,x,graph=NULL)$dens
                          t4     <- plot(self$mdfit.new,x,graph=NULL)$dens
                          res$res <- list(dens=list(t3,t4),results=list(t1,t2))
                        }
                      } else
                      {
                        if (graph == "all") graph = c("doe","densities")
                        if (self$p != 1)
                        {
                          if ("doe" %in% graph)
                          {
                            n        <- combinations(n=self$p,r=2,repeats=FALSE)
                            p        <- list()
                            for (i in 1:nrow(n))
                            {
                              label1 <- n[i,1]
                              label2 <- n[i,2]
                              ind1   <- n[i,1]+self$md$d
                              ind2   <- n[i,2]+self$md$d
                              df     <- data.frame(x=self$doe.init[,ind1],y=self$doe.init[,ind2],col="Initial DOE")
                              df2    <- data.frame(x=self$doe.new[-c(1:nrow(self$doe.init)),ind1],
                                                   y=self$doe.new[-c(1:nrow(self$doe.init)),ind2],col="Additional points")
                              df     <- rbind(df,df2)
                              p[[i]] <- ggplot(df, aes(x = x,y = y, color=col)) + geom_point(shape=3) + theme_light() +
                                xlab(substitute(theta[label1])) + ylab(substitute(theta[label2])) +
                                theme(legend.title = element_blank(),legend.position = c(0.2,0.8),
                                      legend.background = element_rect(linetype="solid", colour ="grey"))
                            }
                            grid.arrange2 <- function(...) return(grid.arrange(...,nrow=1))
                            do.call(grid.arrange2,p)
                            res$doe <- p
                          }
                          if ("densities" %in% graph)
                          {
                            grid.arrange2 <- function(...) return(grid.arrange(...,nrow=1))
                            t1     <- plot(self$mdfit,x,graph=NULL)$out + ggtitle("Before sequential design")
                            t2     <- plot(self$mdfit.new,x,graph=NULL)$out + ggtitle("After sequential design")
                            t3     <- plot(self$mdfit,x,graph=NULL)$dens
                            t4     <- plot(self$mdfit.new,x,graph=NULL)$dens
                            arrangeGrob2 <- function(...){arrangeGrob(...,nrow=1,top="Before sequential design")}
                            arrangeGrob3 <- function(...){arrangeGrob(...,nrow=1,top="After sequential design")}
                            grid.arrange(do.call(arrangeGrob2,t3),do.call(arrangeGrob3,t4),ncol=1)
                            res$res <- list(dens=list(t3,t4),results=list(t1,t2))
                          }
                        } else
                        {
                          if ("densities" %in% graph)
                          {
                            grid.arrange2 <- function(...) return(grid.arrange(...,nrow=1))
                            t1       <- plot(self$mdfit,x,graph=NULL)$out + ggtitle("Before sequential design")
                            t2       <- plot(self$mdfit.new,x,graph=NULL)$out + ggtitle("After sequential design")
                            t3     <- plot(self$mdfit,x,graph=NULL)$dens
                            t4     <- plot(self$mdfit.new,x,graph=NULL)$dens
                            arrangeGrob2 <- function(...){arrangeGrob(...,nrow=1,top="Before sequential design")}
                            arrangeGrob3 <- function(...){arrangeGrob(...,nrow=1,top="After sequential design")}
                            grid.arrange(do.call(arrangeGrob2,t3),do.call(arrangeGrob3,t4),ncol=1)
                            res$res <- list(dens=list(t3,t4),results=list(t1,t2))
                          }
                        }
                      }
                      invisible(res)
                    }
)


