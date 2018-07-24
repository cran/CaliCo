#' Generates \code{\link{model.class}} objects.
#'
#' \code{model} is a function that generates a calibration model and the associated likelihood.
#'
#' @details The different statistical models are: \itemize{\item{Model1:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model2:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model3:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' \item{Model4:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' }
#' where \eqn{for i in [1,\dots,n] \epsilon(x_i)~N(0,\sigma^2)}, \eqn{F(.,.)~PG(m_1(.,.),c_1{(.,.),(.,.)})}
#'  and \eqn{\delta(.)~PG(m_2(.),c_2(.,.))}.
#' There is four kind of models in calibration. They are properly defined in [1].
#'
#' To establish a Gaussian process three options are available:
#' \itemize{
#' \item \strong{opt.gp} is an option list containing the parameters to establish the surrogate (only for model2 and model4).
#' \itemize{
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{DOE}{ design of experiments for the surrogate (default value NULL). If NULL the DOE is automatically
#' generated with the \strong{opt.emul} option.}}
#' }
#' \item \strong{opt.emul} is an option list containing characteristics about emulation option (only for model2 and model4).
#' \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (default value 1)}
#' \item{\strong{n.emul}}{ the number of points for contituing the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{binf}{ the lower bound of the parameter vector (default value 0)}}
#' \item{\strong{bsup}{ the upper bound of the parameter vector (default value 1)}}}
#' \item \strong{opt.sim} is an option list containing the design and corresponding outputs of the code, in the case
#' where no numerical code is available (only for model2 and model4).\itemize{
#' \item{\strong{Ysim}}{ Output of the code}
#' \item{\strong{DOEsim}}{ DOE corresponding to the output of the code}}
#' }
#' To add a discrepancy in the model, the option opt.disc must be added:
#' \itemize{
#' \item \strong{opt.disc} is an option list containing characteristics on the discrepancy (only for model3 and model4)
#' \itemize{
#' \item{\strong{kernel.type}{ see \code{\link{kernel.fun}} for further details}}
#' }
#' }
#'
#' @param code the computational code (function of X and theta)
#' @param X the matrix of the forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2","model3","model4")
#' by default "model1" is chosen. See details for further clarifications
#' @param ... additional options (see details section)
#' @return \code{model} returns a \code{model.class} object. This class contains two main methods:
#' \itemize{
#' \item{plot(model,x)}{ this method generates the plot for a new \eqn{\Theta}, \eqn{\Theta_D} (for model3 and model4),
#'  \eqn{\sigma^2}. The parameter values need to be added to the model with the pipe \code{\%<\%}. The argument \code{x}
#'  represents the x-axis and have to be specified to produce a plot.}
#' \item{print(model)}{ this method presents several pieces of information about the model.}
#' }
#' @author M. Carmassi
#' @seealso \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{forecast}}, \code{\link{sequentialDesign}}
#' @references [1] CARMASSI, Mathieu, BARBILLON, Pierre, KELLER, Merlin, et al. Bayesian calibration of
#'  a numerical code for prediction. arXiv preprint arXiv:1801.01810, 2018.
#' @examples
#' \dontrun{
#' ###### The code to calibrate
#' X <- cbind(seq(0,1,length.out=5),seq(0,1,length.out=5))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(5,0,0.1)
#'
#' ###### For the first model
#' ### Generate the model
#' model1 <- model(code,X,Yexp,"model1")
#' ### Plot the results with the first column of X
#' model1 %<% list(theta=c(1,1,11),var=0.01)
#' plot(model1,X[,1],CI="err")
#'
#' ### Summury of the foo class generated
#' print(model1)
#'
#' ###### For the second model
#' ### code function is available, no DOE generated upstream
#' binf <- c(0.9,0.9,10.5)
#' bsup <- c(1.1,1.1,11.5)
#' opt.gp <- list(type="matern5_2", DOE=NULL)
#' opt.emul <- list(p=3,n.emul=150,binf=binf,bsup=bsup,type="maximinLHS")
#' model2 <- model(code,X,Yexp,"model2",
#'                 opt.gp=opt.gp,
#'                 opt.emul=opt.emul)
#' model2 %<% list(theta=c(1,1,11),var=0.1)
#' ### Plot the model
#' plot(model2,X[,1])
#'
#' ### code function is available and use a specific DOE
#' DOE <- DiceDesign::lhsDesign(20,5)$design
#' DOE[,3:5] <- unscale(DOE[,3:5],binf,bsup)
#'
#' opt.gp <- list(type="matern5_2", DOE=DOE)
#' model2 <- model(code,X,Yexp,"model2",
#'                 opt.gp=opt.gp)
#' model2 %<% list(theta=c(1,1,11),var=0.1)
#' plot(model2,X[,1])
#'
#' ### no code function but DOE and code output available
#' Ysim <- matrix(nr=20,nc=1)
#' for (i in 1:20)
#' {
#'   covariates <- as.matrix(DOE[i,1:2])
#'   dim(covariates) <- c(1,2)
#'   Ysim[i] <- code(covariates,DOE[i,3:5])
#' }
#'
#' opt.sim <- list(Ysim=Ysim,DOEsim=DOE)
#' opt.gp <- list(type="matern5_2", DOE=NULL)
#' model2 <- model(code=NULL,X,Yexp,"model2",
#'                 opt.gp=opt.gp,
#'                 opt.sim=opt.sim)
#' model2 %<% list(theta=c(1,1,11),var=0.1)
#' plot(model2,X[,1])
#'
#' ###### For the third model
#' model3 <- model(code,X,Yexp,"model3",opt.disc=list(kernel.type="gauss"))
#' model3 %<% list(theta=c(1,1,11),thetaD=c(20,0.5),var=0.1)
#' plot(model3,X[,1],CI="err")
#' print(model3)
#'
#'}
#'
#'@export
model <- function(code,X,Yexp,model="model1",...)
{
  switch(model,
         model1={
           obj = model1.class$new(code,X,Yexp,model)
           return(obj)
         },
         model2={
           obj = model2.class$new(code,X,Yexp,model,...)
           return(obj)
         },
         model3={
           obj = model3.class$new(code,X,Yexp,model,...)
           return(obj)
         },
         model4={
           obj = model4.class$new(code,X,Yexp,model,...)
           return(obj)
         }
  )
}


#' Generates \code{\link{prior.class}} objects.
#'
#' \code{prior} is a function that generates a \code{\link{prior.class}} containing information about one or
#' several priors. When several priors are selected, the function \code{prior}
#'  returns a list of \code{\link{prior.class}}.
#'
#' @details The densities implemented are defined as follow
#' \itemize{
#' \item{The Gaussian density:
#' \deqn{f(x)=1/(\sigma*\sqrt(2\pi))exp{-1/2*((x-\mu)/\sigma)^2}}
#' where \strong{\eqn{\mu}} and \strong{\eqn{\sigma}} (the mean and the standard deviation)
#' are the two hyperparameters. The vector \eqn{c(\mu,\sigma^2)} is the one looked for in opt.prior.}
#' \item{The Gamma density:
#' \deqn{f(x)=1/(k^a*\Gamma(a))*x^(a-1)*exp(-(x/k))}
#' where \strong{\eqn{a}} and \strong{\eqn{k}} (the shape and the scale)
#' are the two hyperparameters. The vector \eqn{c(a,k)} is the one looked for in opt.prior.}
#' \item{The Uniform density:
#' \deqn{f(x)=1/(b-a)}
#' where \strong{\eqn{a}} and \strong{\eqn{b}} (the upper and the lower bound)
#' are the two hyperparameters. The vector \eqn{c(a,b)} is the one looked for in opt.prior.}
#' }
#'
#' @param  type.prior the vector of the prior types selected. For example type.prior=c("gaussian","gamma")
#' @param opt.prior list of the hyperparameters relatives to the prior selected. If the first prior selected is
#' Gaussian, the hyperparameters would be the mean and the standard deviation. See Details for further clarifications.
#' @return \code{prior} returns a \code{\link{prior.class}} object. Two main methods are available:
#' \itemize{\item{plot()}{ display the probability density of the prior}
#' \item{print()}{ returns the main information concerning the prior distribution}}
#' @author M. Carmassi
#' @seealso \code{\link{model}}, \code{\link{calibrate}}, \code{\link{forecast}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' #### Only one prior is wanted
#' ## For a Gaussian Prior
#' gaussian <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
#' plot(gaussian)
#'
#' #### For several priors
#' priors <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
#' plot(priors$Prior1)
#' plot(priors$Prior2)
#'}
#'@export
prior <- function(type.prior,opt.prior)
{
  n <- length(type.prior)
  if (n == 1)
  {
    switch(type.prior,
           gaussian = {
             obj = gaussian.class$new(type.prior,opt.prior)
             return(obj)
           },
           gamma={
             obj = gamma.class$new(type.prior,opt.prior)
             return(obj)
           },
           unif={
             obj = unif.class$new(type.prior,opt.prior)
             return(obj)
           }
    )
  } else
  {
    NAmes <- c("Prior1")
    res <- list()
    for (i in 1:n)
    {
      if (i>1){NAmes <- cbind(NAmes,paste("Prior",i,sep=""))}
      switch(type.prior[i],
             gaussian = {
               obj = gaussian.class$new(type.prior[i],opt.prior[[i]])
             },
             gamma={
               obj = gamma.class$new(type.prior[i],opt.prior[[i]])
             },
             unif={
               obj = unif.class$new(type.prior[i],opt.prior[[i]])
             }
      )
      res[[i]] <- obj
      names(res) <- NAmes
    }
    return(res)
  }
}


#' Generates \code{\link{calibrate.class}} objects
#'
#' \code{calibrate} is a function that allows to generate a \code{\link{calibrate.class}}
#'  class in which the estimation is
#' done for a defined \code{\link{model.class}} and \code{\link{prior.class}} objects.
#'
#' @useDynLib CaliCo
#'
#' @param md a \code{\link{model.class}} object
#' @param pr a \code{\link{prior.class}} object
#' @param opt.estim estimation options \itemize{\item{Ngibbs}{ Number of iteration of the algorithm
#' Metropolis within Gibbs}
#' \item{Nmh}{ Number of iteration of the Metropolis Hastings algorithm}
#' \item{thetaInit}{ Initial point}
#' \item{r}{ regulation percentage in the modification of the k in the Metropolis Hastings}
#' \item{sig}{ Covariance matrix for the proposition distribution (\eqn{k*sig})}
#' \item{Nchains}{ Number of MCMC chains to run (if Nchain>1 an output is created called mcmc which
#'  is a coda object \code{\link{codamenu}})}
#' \item{burnIn}{ Number of iteration to withdraw}
#' }
#' @param opt.valid list of cross validation options (default value opt.valid=NULL)\itemize{
#' \item{nCV}{ Number of iterations for the cross validation}
#' \item{type.valid}{ Type of cross validation selected. "loo" (leave one out) is the only method
#'  emplemented so far.}
#' }
#' @return \code{calibrate} returns a \code{\link{calibrate.class}} object. Two main methods are available:
#' \itemize{\item{plot(mdfit, x, graph)}{ displays a series of graphs (ACF, MCMC, density a priori vs a posteriori
#' , correlation between parameters, results on the quantify of interest, etc..) or return a list with all
#' the graphs:}
#' \itemize{
#' \item {mdfit}{ The calibrated model (a \code{\link{calibrate.class}} object)}
#' \item {x}{ The x-axis}
#' \item {graph}{ Allows to select the wanted display. By default all the layout pannel graphs are displayed and
#' \code{graph="all"}. If \code{graph="chains"}, only the layout of the autocorrelation, chains points and densities
#'  a priori and a posteriori is produced. If \code{graph="corr"}, only the layout of the correlation graph between
#'   each parameter is displayed. If \code{graph="result"}, only the result on the quantity of interest is given.
#'   If \code{graph=NULL}, no graphs are produced automatically.}}
#' \item{print(mdfit)}{ returns the main information concerning the \code{\link{calibrate.class}} object}}
#' @author M. Carmassi
#' @seealso \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{forecast}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' ###################### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#' ############### For the first model
#' ###### Definition of the model
#' md <- model(code,X,Yexp,"model1")
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=200,Nmh=400,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=1,burnIn=100)
#' ###### Run the calibration
#' mdfit <- calibrate(md,pr,opt.estim)
#' ####### The plot generated is a list of ggplot
#' p <- plot(mdfit,X[,1])
#' p$out
#' print(mdfit)
#'}
#' @export
calibrate <-function(md,pr,opt.estim,opt.valid=NULL)
{
  res <- calibrate.class$new(md,pr,opt.estim,opt.valid,onlyCV=FALSE)
  return(res)
}

#' Generates a forecast base on calibration run with \code{\link{calibrate}}
#'
#' \code{forecast} is a function that allows to generate a new \code{\link{model.class}} in which the prediction is
#' done with the Maximum A Posteriori
#'
#' Note that all the methods for a \code{\link{model.class}} object are availble. Be careful with the \code{x} in the
#'  plot function. It needs to be the x-axis of calibrated data and predicted data.
#'
#' @useDynLib CaliCo
#'
#' @param modelfit a \code{\link{calibrate.class}} object
#' @param x.new newdata for the prediction
#' @return return a \code{\link{model.class}} (see \code{\link{model.class}} for more details)
#' @author M. Carmassi
#' @seealso \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' ###################### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#' ############### For the first model
#' ###### Definition of the model
#' md <- model(code,X,Yexp,"model1")
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=200,Nmh=600,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=1,burnIn=100)
#' ###### Run the calibration
#' mdfit <- calibrate(md,pr,opt.estim)
#' ###### Prediction between 1 and 1.2
#' X.new <- cbind(seq(1,1.2,length.out=10),seq(1,1.2,length.out=10))
#' fr <- forecast(mdfit,X.new)
#' print(fr)
#' plot(fr,c(X[,1],X.new[,1]))
#'}
#' @export
forecast <-function(modelfit,x.new)
{
  if ("seqDesign.class" %in% class(modelfit))
  {
   md <- modelfit$md.new
   modelfit <- modelfit$mdfit.new
  } else md <- modelfit$md

  if (is.matrix(x.new)) X <- rbind(md$X,x.new)
  else X <- c(md$X,x.new)
  if (md$model %in% c("model1", "model3"))
  {
    md.new <- model(code = md$code, X = X, Yexp = md$Yexp, model = md$model,opt.disc=md$opt.disc)
  } else
  {
    md.new <- model(code = md$code, X = X, Yexp = md$Yexp, model = md$model,opt.disc=md$opt.dis, opt.gp= md$opt.gp,
                    opt.emul = md$opt.emul, opt.sim = md$opt.sim)
  }
  l <- length(modelfit$output$MAP)
  if (md$model %in% c("model1","model2"))
  {
    options(warn = -1)
    md.new %<% list(theta=modelfit$output$MAP[-l],var=modelfit$output$MAP[l])
    options(warn = 0)
  } else
  {
    options(warn = -1)
    md.new %<% list(theta=modelfit$output$MAP[-((l-2):l)],thetaD=modelfit$output$MAP[((l-2):(l-1))],
                    var=modelfit$output$MAP[l])
    options(warn = 0)
  }
  fr <- forecast.class$new(modelfit,md.new,x.new)
  return(fr)
}


#' Calibration with a sequential design
#'
#' The aim is to reduce the
#' error produced by the initial estimation of the Gaussian process by fortifying the initial DOE. The method consists
#' in proposing new points based on the expectancy improvement criterion. The method and the algorithm are detailed in
#' [Damblin et al. 2018]
#'
#' @param md the model to improve (model2 or model4)
#' @param pr list of priors to use for calibration
#' @param opt.estim estimation options \itemize{\item{Ngibbs}{Number of iteration of the algorithm
#' Metropolis within Gibbs}
#' \item{Nmh}{ Number of iteration of the Metropolis Hastings algorithm}
#' \item{thetaInit}{ Initial point}
#' \item{r}{ regulation percentage in the modification of the k in the Metropolis Hastings}
#' \item{sig}{ Covariance matrix for the proposition distribution (\eqn{k*sig})}
#' \item{Nchains}{ Number of MCMC chains to run (if Nchain>1 an output is created called mcmc which
#'  is a coda object \code{\link{codamenu}})}
#' \item{burnIn}{ Number of iteration to withdraw}
#' }
#' @param k number of iteration in the algorithm
#' @return a \code{\link{seqDesign.class}}
#' @seealso \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' ###### The code to calibrate
#' X <- cbind(seq(0,1,length.out=5),seq(0,1,length.out=5))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(5,0,0.1)
#'
#' ###### For the second model
#' ### code function is available, no DOE generated upstream
#' binf <- c(0.9,0.9,10.5)
#' bsup <- c(1.1,1.1,11.5)
#' opt.gp <- list(type="matern5_2", DOE=NULL)
#' opt.emul <- list(p=3,n.emul=150,binf=binf,bsup=bsup,type="maximinLHS")
#' model2 <- model(code,X,Yexp,"model2",
#'                 opt.gp=opt.gp,
#'                 opt.emul=opt.emul)
#' model2 %<% list(theta=c(1,1,11),var=0.1)
#'
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=200,Nmh=400,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=1,burnIn=100)
#' ###### Run the sequential calibration
#' mdfit <- sequentialDesign(model2,pr,opt.estim,2)
#' #plot(mdfit,X[,1])
#' }
#' @references DAMBLIN, Guillaume, BARBILLON, Pierre, KELLER, Merlin, et al. Adaptive numerical designs for the
#'  calibration of computer codes. SIAM/ASA Journal on Uncertainty Quantification, 2018, vol. 6, no 1, p. 151-179.
#' @export
sequentialDesign <- function(md,pr,opt.estim,k)
{
  obj <- seqDesign.class$new(md,pr,opt.estim,k)
  return(obj)
}



#' Return Maximum A Posteriori (MAP) and Mean A Posteriori estimation of a calibration
#'
#' \code{estimators} is a function that returns a list of two elements which are the MAP and the Mean A Posteriori
#'
#' @useDynLib CaliCo
#'
#' @param modelfit a \code{\link{calibrate.class}} object
#' @return return a \code{list}:
#' \itemize{
#' \item {MAP}{ The Maximum A Posteriori}
#' \item {MEAN}{ The Mean A Posteriori}}
#' @author M. Carmassi
#' @seealso \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' ###################### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#' ############### For the first model
#' ###### Definition of the model
#' md <- model(code,X,Yexp,"model1")
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=200,Nmh=600,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=1,burnIn=100)
#' ###### Run the calibration
#' mdfit <- calibrate(md,pr,opt.estim)
#' estimators(mdfit)
#' }
#' @export
estimators <-function(modelfit)
{
  MEAN <- apply(modelfit$output$out$THETA[-c(1:modelfit$opt.estim$burnIn),],2,mean)
  return(list(MAP=modelfit$output$MAP,MEAN=MEAN))
}


#' Return the MCMC chain in a \code{data.frame}
#'
#' \code{chain} is a function that returns a \code{data.frame} of calibration run without the burn-in
#'
#' @useDynLib CaliCo
#'
#' @param modelfit a \code{\link{calibrate.class}} object
#' @param coda if TRUE returns a coda object (if Nchains in opt.estim is higher than 1 a coda object
#'  is automatically returned see \code{\link{codamenu}})
#' @return return a \code{data.frame} or a coda object of the MCMC chain(s) generated.
#' @author M. Carmassi
#' @seealso \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}}, \code{\link{sequentialDesign}}
#' @examples
#' \dontrun{
#' ###################### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#' ############### For the first model
#' ###### Definition of the model
#' md <- model(code,X,Yexp,"model1")
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=200,Nmh=600,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=1,burnIn=100)
#' ###### Run the calibration
#' mdfit <- calibrate(md,pr,opt.estim)
#'
#' mcmc <- chain(mdfit)
#' ### Check coda object
#' is.mcmc(mcmc)
#' ### get all the chain
#' mcmc <- chain(mdfit,coda=FALSE)
#' head(mcmc)
#'
#' ### Multi chains
#' opt.estim=list(Ngibbs=200,Nmh=600,thetaInit=c(1,1,11,0.1),r=c(0.3,0.3),
#' sig=diag(4),Nchains=2,burnIn=100)
#' ###### Run the calibration
#' mdfit2 <- calibrate(md,pr,opt.estim)
#'
#' mcmc <- chain(mdfit2)
#' is.mcmc.list(mcmc)
#' }
#' @export
chain <-function(modelfit,coda=TRUE)
{
  if (modelfit$opt.estim$Nchain == 1)
  {
    if (coda == TRUE) return(modelfit$mcmc)
    else
    {
      cc <- modelfit$output$out$THETA[-c(1:modelfit$opt.estim$burnIn),]
      L <- ncol(cc)
      if (modelfit$md$model %in% c("model1", "model2")){
        p <- L-1
        Names <- NULL
        for (i in 1:p) Names <- c(Names,paste("theta",i,sep="_"))
        Names <- c(Names,"Var")
      } else {
        p <- L-3
        for (i in 1:p) Names <- c(Names,paste("theta",i,sep="_"))
        Names <- c(Names,"sigD","psiD","Var")
      }
      colnames(cc) <- Names
      return(as.data.frame(cc))
    }
  } else
  {
    return(modelfit$mcmc)
  }
}


#' Generates covariances matrices thanks to \code{\link{Kernel.class}}
#'
#' \code{Kernel.fun} is a function that allows us to generate covariances matrices from data
#'
#' @param X data
#' @param var the variance for the covariance function
#' @param psi the parameter vector
#' @param kernel.type the choice of the form of the kernel (with d chosen as an euclidian distance) \itemize{
#' \item{gauss}{ \deqn{\sigma^2 exp{-1/2(d/\psi)^2}}}
#' \item{exp}{ \deqn{\sigma^2 exp{-1/2 d/\psi}}}
#' \item{matern3_2}{ \deqn{\sigma^2(1+\sqrt{3}d^2/\psi) exp{-\sqrt{3}d^2/\psi}}}
#' \item{matern5_2}{ \deqn{\sigma^2(1+\sqrt{5}d^2/\psi+5d^2/(3\psi^2))exp{-\sqrt{5}d^2/\psi}}}}
#' @return \code{Kernel.fun} returns a covariance matrix
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{prior.class}}
#' @examples
#' \dontrun{
#' X <- cbind(seq(0,10,length.out=10),seq(8,20,length.out=10))
#' var <- 2
#' psi <- 0.1
#' Cov <- kernel.fun(X,var,psi,kernel.type="matern5_2")
#'}
#' @export
kernel.fun <- function(X,var,psi,kernel.type="gauss")
{
  if (is.null(kernel.type)){kernel.type <- "gauss"}
  switch(kernel.type,
         gauss={
           obj = gauss.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         exp={
           obj = exp.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         matern3_2={
           obj = matern3_2.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         matern5_2={
           obj = matern5_2.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         }
  )
}

#' Function which unscale a vector between two bounds
#'
#' @param  x the vector to unscale
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @return y the vector unscaled
#' @examples
#' \dontrun{
#' X <- runif(3)
#' Y <-unscale.vector(X,rep(10,3),rep(15,3))
#' print(Y)
#' }
#' @export
unscale.vector <- function(x,binf,bsup){
  n = length(x)
  if (length(binf)!=length(bsup))
  {
    print("Please enter bounds of the same size")
  }
  if (length(binf)==1)
  {
    y=rep(1,n)
    for (i in 1:n){
      y[i]<- bsup*x[i]+(1-x[i])*binf
    }
  }else
  {
    y=rep(1,n)
    for (i in 1:n){
      y[i]<- bsup[i]*x[i]+(1-x[i])*binf[i]
    }
  }
  return(y)
}


#' Function which unscale only the diagonal component of a matrix
#'
#' @param  M the matrix
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @return the normalized diagonal
#' @examples
#' \dontrun{
#' X <- diag(3)*runif(3)
#' Y <-unscale.matrix.diag(X,rep(10,3),rep(15,3))
#' print(Y)
#' }
#' @export
unscale.matrix.diag <- function(M,binf,bsup){
  n <- dim(M)[1]
  T <- rep(1,n)
  for (i in 1:n){
    T[i] <- M[i,i]
  }
  T <- unscale.vector(T,binf,bsup)
  for (i in 1:n){
    M[i,i] <- T[i]
  }
  return(M)
}


#' Function which unscale un matrix or a vector
#'
#' @param M the matrix or the vector
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @param diag default value False if we want to unscale the whole matrix
#' @param sym default value False if we do not have a symmetric matrix
#' @return the unscaled vector or matrix
#' @examples
#' \dontrun{
#' X <- diag(3)*runif(3)
#' Y <- unscale(X,rep(10,3),rep(15,3))
#' print(Y)
#' }
#' @export
unscale <- function(M,binf,bsup,diag=FALSE,sym=FALSE){
  if (diag==FALSE){
    if(sym==FALSE){
      if (is.matrix(M)==FALSE){
        if (length(binf)==length(M) & length(bsup)==length(M)){
          temp <- unscale.vector(M,binf,bsup)
        } else {
          temp <- unscale.vector(M,binf,bsup)
        }
      } else{
        N <- nrow(M)
        n <- ncol(M)
        if (n==1){
          temp <- unscale.vector(M,binf,bsup)
        } else {
          if (n!=length(binf) & N!=length(bsup)){
            print('The bound s sizes are not compatible.
                  Please enter a upper bound and a lower bound of the same size than the matrix.')
          } else {
            temp <- matrix(nrow = N,ncol = n)
            for (i in 1:N){
              temp[i,] <- unscale.vector(M[i,],binf,bsup)
            }
          }
        }
    }
    } else {
      # Pour les matrice symétriques en cours...
    }
    } else {
      temp<-unscale.matrix.diag(M,binf,bsup)
    }
  return(temp)
  }


#' Function that deals with negative eigen values in a matrix not positive definite
#'
#' @param X the matrix or the vector
#' @return the new matrix X
#' @export
DefPos <- function(X)
{
  p <- eigen(X)$vectors
  e <- eigen(X)$values
  if (all(e>0)){} else
  {
    e[which(e<0)] <- 1e-4
  }
  d <- diag(e)
  return(t(p)%*%d%*%p)
}

#' Simulate from a Multivariate Normal Distribution
#'
#' The matrix decomposition is done via eigen; although a Choleski decomposition might be faster,
#' the eigen decomposition is stabler.
#'
#' @param n the number of samples required.
#' @param mu a vector giving the means of the variables.
#' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix of the variables.
#' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
#' @param empirical logical. If true, mu and Sigma specify the empirical not population mean and covariance matrix.
#' @param EISPACK logical: values other than FALSE are an error.
#' @return If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.
#' @export
multivariate <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}


#' Operator to define active bindings variables
#' @param md a statistical model defined by the function \code{\link{model}}
#' @param param a \code{list} of parameter values
#' @return a \code{\link{model.class}} parametrized.
#' @export
"%<%" <- function (md,param)
{
  if ("model.class" %in% class(md))
  {
    if ("theta" %in% names(param) & "var" %in% names(param))
    {
      warning("Please be carefull to the size of the parameter vector",call. = FALSE)
      if (length(param$var) > 1)
      {
        stop("Wrong variance size",call. = FALSE)
      }
      md$theta    <- param$theta
      md$var      <- param$var
    } else
    {
      stop("To realize a parametrization of the model please enter a list containing theta and var",call. = FALSE)
    }
    if (md$model %in% c("model3","model4"))
    {
      if ("theta" %in% names(param) & "var" %in% names(param))
      {
        if (length(param$thetaD) != 2)
        {
          stop("Wrong discrepancy parameter size",call. = FALSE)
        }
        md$thetaD <- param$thetaD
      } else
      {
        stop("For the third model 3 and 4, thetaD has to be added",call. = FALSE)
      }
    }
  } else
  {
    stop("Not a model.class",call. = FALSE)
  }
}





