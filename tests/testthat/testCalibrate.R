x <- seq(0,1,length.out = 5)
code <- function(x,theta){return(x*theta)}
Yexp <- code(x,1)
md1 <- model(code = code, X=x,Yexp = Yexp,model = "model1")
binf <- 0.5
bsup <- 1.5
opt.emul <- list(p=1,n.emul=10,type="matern5_2",binf=binf,bsup=bsup,DOE=NULL)
md2 <- model(code = code, X=x,Yexp = Yexp,model = "model2",opt.emul = opt.emul)
opt.disc=list(kernel.type="matern5_2")
md3 <- model(code = code, X=x,Yexp = Yexp,model = "model3",opt.disc = opt.disc)
md4 <- model(code = code, X=x,Yexp = Yexp,model = "model3",opt.disc = opt.disc,opt.emul = opt.emul)
pr1 <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(1,0.01),c(1,0.01)))
pr2 <- prior(type.prior=c("gaussian","unif","gamma","gamma"),opt.prior=list(c(1,0.01),c(0,1),c(2,0.1),c(1,0.01)))

#calibrate testing
opt.estim1=list(Ngibbs=200,Nmh=600,thetaInit=c(1,0.01),k=c(6e-3,1e-4),sig=diag(2),Nchains=1,burnIn=300)
opt.estim2=list(Ngibbs=200,Nmh=600,thetaInit=c(1,0.5,0.01,0.01),k=rep(5e-3,4),sig=diag(4),Nchains=1,burnIn=300)
mdfit1 <- calibrate(md1,pr1,opt.estim1)
mdfit2 <- calibrate(md2,pr1,opt.estim1)
mdfit3 <- calibrate(md3,pr2,opt.estim2)
mdfit4 <- calibrate(md4,pr2,opt.estim2)
expect_is(mdfit1,"calibrate.class")
expect_is(mdfit2,"calibrate.class")
expect_is(mdfit3,"calibrate.class")
expect_is(mdfit4,"calibrate.class")
