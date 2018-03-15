## ---- fig.width=4,echo=FALSE,fig.align='center'--------------------------
library(ggplot2)

n=10
X <- seq(0,1,length.out=n)
code <- function(X,theta)
{
  return((6*X-2)^2*sin(theta*X-4))
}
Yexp <- code(X,10.5)+rnorm(n,0,sqrt(0.01))
Yt <- code(X,11)
gdata <- data.frame(y=Yexp,x=X,type="Experiments")
gdata2 <- data.frame(y=Yt,x=X,type="Theoretical results")
gdata <- rbind(gdata,gdata2)

p1 <- ggplot(gdata,aes(x = x,y=y, col=type))+geom_line() + theme_light() + theme(legend.position=c(0.50,0.80),
                                                                                 legend.title=element_blank())
p1

## ---- fig.width=6,echo=FALSE,fig.align='center'--------------------------
theta.prior <- rnorm(100,11,sqrt(.5))
sigma.prior <- 0.01

Y_prior <- matrix(nr=100,nc=n)
for (i in 1:100){
  Y_prior[i,] <- code(X,theta.prior[i])+rnorm(n,0,sqrt(sigma.prior))
}

qq <- apply(Y_prior,2,quantile,probs=c(0.05,0.95))
gdata <- data.frame(y=Yexp,x=X,upper=qq[2,],lower=qq[1,],type="Experiments",fill="credibility interval a priori at 90%")
gdata2 <- data.frame(y=Yt,x=X,upper=qq[2,],lower=qq[1,],type="Prior belief",fill="credibility interval a priori at 90%")
gdata <- rbind(gdata,gdata2)

p2 <- ggplot(gdata)+geom_line(aes(x=x,y=y,col=type))+geom_ribbon(aes(x=x,ymax=upper,ymin=lower,fill=fill),alpha=0.3)+theme_light()+theme(legend.title=element_blank(),legend.key=element_rect(colour=NA),legend.text=element_text(size = '8'))
p2

## ---- echo=TRUE----------------------------------------------------------
library(CaliCo)
# Number of experiments
n=10
# Time interval
t <- seq(0,1,length.out=n)
# Code definition
code <- function(t,theta)
{
  return((6*t-2)^2*sin(theta*t-4))
}
# Generate the experiment
Yexp <- code(t,10.5)+rnorm(n,0,sqrt(0.01))
# Generate the first model
model1 <- model(code=code,X=t,Yexp=Yexp,model="model1")

## ----echo=TRUE-----------------------------------------------------------
print(model1)

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
plot(model1,11,0.01)

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
library(ggplot2)
p <- plot(model1,11,0.01)
p+ggtitle("Model1 and experiments")+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ----echo=TRUE-----------------------------------------------------------
# Set the lower and upper bound (taken large intentionally)
binf <- 8
bsup <- 13

# Set the emulation option
opt.emul <- list(p=1,n.emul=40,type="matern3_2",binf=binf,bsup=bsup,DOE=NULL)
# Generate the second model
model2 <- model(code,X,Yexp,"model2",opt.emul)

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
p <- plot(model2,11,0.01)
p+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ---- echo=TRUE----------------------------------------------------------
opt.disc <- list(kernel.type="matern5_2")
# Generate the third model
model3 <- model(code,X,Yexp,"model3",opt.disc=opt.disc)

# Generate the fourth model where the emulation option are needed
### Here opt.emul is the same as the one for the second model
model4 <- model(code,X,Yexp,"model4",opt.emul,opt.disc)

## ---- echo=TRUE, fig.width=4, fig.align="center"-------------------------
p <- plot(model3,11,c(2,0.5),0.01)
p+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ---- echo=TRUE, fig.width=6, fig.align="center"-------------------------
p <- plot(model4,11,c(1,0.1),0.01)
p+ylab("y")+xlab("t")+
  theme(legend.position="right",legend.text=element_text(size = '7'))

## ---- fig.width=4,fig.align='center'-------------------------------------
gaussian <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
plot(gaussian)

## ---- fig.width=4,fig.align='center'-------------------------------------
priors <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
plot(priors$Prior1)
plot(priors$Prior2)

## ---- echo=TRUE----------------------------------------------------------
print(priors$Prior1)

## ---- echo=TRUE----------------------------------------------------------
print(priors$Prior2)

## ---- echo=TRUE----------------------------------------------------------
# Definition of the estimation options for the model 1 and model 2
opt.estim1=list(Ngibbs=1000,Nmh=3000,thetaInit=c(11,0.01),k=c(1e-4,1e-4),sig=diag(2),Nchains=1,burnIn=1000)
opt.estim2=list(Ngibbs=2000,Nmh=4000,thetaInit=c(11,0.01),k=c(1e-5,1e-5),sig=diag(2),Nchains=1,burnIn=1000)

# Generate the prior densities
pr1 <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(11,0.5),c(2,0.01)))

# Run the Bayesian Calibration for the first
mdfit1 <- calibrate(model1,pr1,opt.estim1)

## ---- echo=TRUE----------------------------------------------------------
# Run the Bayesian Calibration for the first
mdfit2 <- calibrate(model2,pr1,opt.estim2)

## ---- echo=TRUE----------------------------------------------------------
print(mdfit1)

## ---- echo=TRUE----------------------------------------------------------
print(mdfit2)

## ---- echo=TRUE, fig.width=4---------------------------------------------
plot(mdfit1)[[4]]

## ---- echo=TRUE, fig.width=6---------------------------------------------
Yt <- code(X,11)
gdata <- data.frame(y=Yt,x=X,type="Theoretical results")

plot(mdfit1)[[4]]+geom_line(aes(x,y,colour=type), gdata)+ggtitle("Results for Model1")+xlab("x")+ylab("y")+theme(legend.position="right",legend.text=element_text(size = '7'))

## ---- echo=TRUE, fig.width=6---------------------------------------------
t2 <- mdfit2$plot()
Yt <- code(X,11)
gdata <- data.frame(y=Yt,x=X,type="Theoretical results")

t2[[4]]+geom_line(aes(x,y,colour=type), gdata)+ggtitle("Results for Model2")+xlab("x")+ylab("y")+theme(legend.position="right",legend.text=element_text(size = '7'))

## ---- echo=TRUE----------------------------------------------------------
# Definition of the estimation options for the model 3 and model 4
opt.estim2=list(Ngibbs=1000,Nmh=3000,thetaInit=c(11,0.5,0.01,0.01),k=c(1e-5,1e-6,1e-7,1e-7),sig=diag(4),Nchains=1,burnIn=1000)

# Generate the prior densities
pr2 <- prior(type.prior=c("gaussian","gamma","unif","gamma"),opt.prior=list(c(11,0.5),c(1,0.5),c(0,1),c(2,0.01)))

# Run the Bayesian Calibration for the third
mdfit3 <- calibrate(model3,pr2,opt.estim2)

## ---- echo=TRUE----------------------------------------------------------
# Run the Bayesian Calibration for the fourth model
mdfit4 <- calibrate(model4,pr2,opt.estim2)

## ---- echo=TRUE, fig.width=6---------------------------------------------
t3 <- plot(mdfit3)
t4 <- plot(mdfit4)
Yt <- code(X,11)
gdata <- data.frame(y=Yt,x=X,type="Theoretical results")

t3[[4]]+geom_line(aes(x,y,colour=type), gdata)+ggtitle("Results for Model3")+xlab("x")+ylab("y")+theme(legend.position="right",legend.text=element_text(size = '7'))
t4[[4]]+geom_line(aes(x,y,colour=type), gdata)+ggtitle("Results for Model4")+xlab("x")+ylab("y")+theme(legend.position="right",legend.text=element_text(size = '7'))

