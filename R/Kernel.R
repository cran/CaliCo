#' A Reference Class to generates differents model objects
#'
#' @description See the function \code{\link{model}} which produces an instance of this class
#' @field code a function which takes in entry X and theta
#' @field X the matrix of the forced variables
#' @field var the variance for the covariance function
#' @field psi the scale vector of correlation length
#' @field n number of experiments
#' @field Kernel.type the chosen form of covariance
#' @field Cov the covariance
#' @export
Kernel.class <- R6Class(classname = "Kernel.class",
                           public = list(
                             X           = NULL,
                             var         = NULL,
                             psi         = NULL,
                             n           = NULL,
                             Kernel.type = NULL,
                             Cov         = NULL,
                             initialize = function(X=NA,var=NA,psi=NA,Kernel.type=NA)
                             {
                               self$X           <- X
                               self$var         <- var
                               self$psi         <- psi
                               if (is.null(dim(self$X))==TRUE)
                               {
                                 self$n <- length(self$X)
                               } else
                               {
                                 self$n <- nrow(self$X)
                               }
                               self$Kernel.type <- Kernel.type
                             }
                           ))


gauss.class <- R6Class(classname= "gauss.class",
                           inherit = Kernel.class,
                           public = list(
                             initialize = function(X,var,psi,Kernel.type="gauss")
                             {
                               super$initialize(X,var,psi,Kernel.type="gauss")
                               self$Cov <- self$covariance()
                             },
                             covariance = function()
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               if (is.null(ncol(self$X)))
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 for (j in 1:self$n)
                                 {
                                   for (i in 1:self$n)
                                   {
                                     self$Cov[i,j] <- self$var*exp(-1/2*(sum((self$X[i]-self$X[j])^2)/self$psi)^2)
                                   }
                                 }
                               } else
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 for (j in 1:self$n)
                                 {
                                   for (i in 1:self$n)
                                   {
                                     self$Cov[i,j] <- self$var*exp(-1/2*(sum((self$X[i,]-self$X[j,])^2)/self$psi)^2)
                                   }
                                 }
                               }
                               return(self$Cov)
                             }
                           ))

exp.class <- R6Class(classname= "exp.class",
                           inherit = Kernel.class,
                         public = list(
                           initialize = function(X,var,psi,Kernel.type="exp")
                           {
                             super$initialize(X,var,psi,Kernel.type="exp")
                             self$Cov <- self$covariance()
                           },
                           covariance = function(X,var,psi)
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             if (is.null(ncol(self$X)))
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               for (j in 1:self$n)
                               {
                                 for (i in 1:self$n)
                                 {
                                   self$Cov[i,j] <- self$var*exp(sum((self$X[i]-self$X[j])^2)/self$psi)
                                 }
                               }
                             } else
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               for (j in 1:self$n)
                               {
                                 for (i in 1:self$n)
                                 {
                                   self$Cov[i,j] <- self$var*exp(sum((self$X[i,]-self$X[j,])^2)/self$psi)
                                 }
                               }
                             }
                             return(self$Cov)
                           }
))


matern3_2.class <- R6Class(classname= "matern3_2.class",
                         inherit = Kernel.class,
                         public = list(
                         initialize = function(X,var,psi,Kernel.type="matern3_2")
                         {
                           super$initialize(X,var,psi,Kernel.type="matern3_2")
                           self$Cov <- self$covariance()
                         },
                         covariance = function(X,var,psi)
                         {
                           self$Cov <- matrix(nr=self$n,nc=self$n)
                           if (is.null(ncol(self$X)))
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             for (j in 1:self$n)
                             {
                               for (i in 1:self$n)
                               {
                                 self$Cov[i,j] <- self$var*(1+sqrt(3)*sum((self$X[i]-self$X[j])^2)/self$psi)*
                                   exp(-sqrt(3)*sum((self$X[i]-self$X[j])^2)/self$psi)
                               }
                             }
                           } else
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             for (j in 1:self$n)
                             {
                               for (i in 1:self$n)
                               {
                                 self$Cov[i,j] <- self$var*(1+sqrt(3)*sum((self$X[i,]-self$X[j,])^2)/self$psi)*
                                   exp(-sqrt(3)*sum((self$X[i,]-self$X[j,])^2)/self$psi)
                               }
                             }
                           }
                           return(self$Cov)
                         }
))


matern5_2.class <- R6Class(classname= "matern5_2.class",
                               inherit = Kernel.class,
                               public = list(
                               initialize = function(X,var,psi,Kernel.type="matern5_2")
                               {
                                 super$initialize(X,var,psi,Kernel.type="matern5_2")
                                 self$Cov <- self$covariance()
                               },
                               covariance = function(X,var,psi)
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 if (is.null(ncol(self$X)))
                                 {
                                   self$Cov <- matrix(nr=self$n,nc=self$n)
                                   for (j in 1:self$n)
                                   {
                                     for (i in 1:self$n)
                                     {
                                       self$Cov[i,j] <- self$var*(1+sqrt(5)*sum((self$X[i]-self$X[j])^2)/self$psi+
                                                               5/3*(sum((self$X[i]-self$X[j])^2)/self$psi)^2)*
                                                               exp(-sqrt(5)*sum((self$X[i]-self$X[j])^2)/self$psi)
                                     }
                                   }
                                 } else
                                 {
                                   self$Cov <- matrix(nr=self$n,nc=self$n)
                                   for (j in 1:self$n)
                                   {
                                     for (i in 1:self$n)
                                     {
                                       self$Cov[i,j] <- self$var*(1+sqrt(5)*sum((self$X[i,]-self$X[j,])^2)/self$psi+
                                                               5/3*(sum((self$X[i,]-self$X[j,])^2)/self$psi)^2)*
                                                               exp(-sqrt(5)*sum((self$X[i,]-self$X[j,])^2)/self$psi)
                                     }
                                   }
                                 }
                                 return(self$Cov)
                               }
))
