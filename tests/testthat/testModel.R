x <- seq(0,1,length.out = 5)
code <- function(x,theta){return(x*theta)}
Yexp <- code(x,1)
#model testing
md1 <- model(code = code, X=x,Yexp = Yexp,model = "model1")
expect_is(md1,"model.class")
binf <- 0.5
bsup <- 1.5
opt.emul <- list(p=1,n.emul=10,type="matern5_2",binf=binf,bsup=bsup,DOE=NULL)
md2 <- model(code = code, X=x,Yexp = Yexp,model = "model2",opt.emul = opt.emul)
expect_is(md2,"model.class")
opt.disc=list(kernel.type="matern5_2")
md3 <- model(code = code, X=x,Yexp = Yexp,model = "model3",opt.disc = opt.disc)
expect_is(md3,"model.class")
md4 <- model(code = code, X=x,Yexp = Yexp,model = "model3",opt.disc = opt.disc,opt.emul = opt.emul)
expect_is(md4,"model.class")


