###############################
### Imbens sensitivity analysis
### Fake data
### Note: the phrase to input in the secureshell window before trying to run this so that we can close the window is:  nohup R --vanilla<Test1Code.R &
### The host: 128.151.41.4
##########################

source("sensitivity2MLE.R")

load("alpha1.RData")

constant <- 1
vars <- with(test.1,na.omit(cbind(constant,X1,X2,X3,W,Y)))
xmatt <- vars[,1:4]
W <- vars[,5]
Y <- vars[,6]


starting <- function(xmat){
  coef.lm <- glm(Y~xmat+W-1)$coefficients  ##starting values for beta
  coef.logit <- glm(W~xmat-1)$coefficients ##starting values for gamma
  coef.lm <- coef.lm[-length(coef.lm)]
  startv <- c(coef.logit,coef.lm,1) ##starting values for gamma,beta, and tau
  startv
}

stvalnt <- starting(xmatt)

aldelval <- cbind(rep(seq(-10,10,.4),51),rep(seq(-10,10,.4),each=51))

ncomb <- nrow(aldelval)

Izero <- imbensbinary(0,0,xmatt,stvalnt)
Outcoef <- c(Izero@Beta,Izero@Tau)
R2yzero <- (t(Outcoef)%*%Izero@VCW%*%Outcoef)/(t(Outcoef)%*%Izero@VCW%*%Outcoef+((pi^2)/3))
R2wzero <- (t(Izero@Gamma)%*%Izero@VC%*%Izero@Gamma)/(t(Izero@Gamma)%*%Izero@VC%*%Izero@Gamma+((pi^2)/3))

rsquare <- function(i){

  Imb <- imbensbinary(aldelval[i,1],aldelval[i,2],xmatt,stvalnt)
  Tau <- Imb@Tau
  Outcf <- c(Imb@Beta,Imb@Tau)
  R2y <- (t(Outcf)%*%Imb@VCW%*%Outcf+(((aldelval[i,2])^2)/4))/(t(Outcf)%*%Imb@VCW%*%Outcf+(((aldelval[i,2])^2)/4)+((pi^2)/3))
  R2w <- (t(Imb@Gamma)%*%Imb@VC%*%Imb@Gamma+(((aldelval[i,1])^2)/4))/(t(Imb@Gamma)%*%Imb@VC%*%Imb@Gamma+(((aldelval[i,1])^2)/4)+((pi^2)/3))
  R2p <- (R2y-R2yzero)/(1-R2yzero)
  R2wp <- (R2w-R2wzero)/(1-R2wzero)
  R2p <- round(R2p,digits=3)
  R2wp <- round(R2wp,digits=3)
  vals <- cbind(Tau, R2p, R2wp, R2y)
  vals
  
}
  
usecluster <- TRUE #if using the cluster enter "TRUE"

if (usecluster) { 
  library(snow)
  nclus <- 24
  clus1 <- makePVMcluster(nclus)
  clusterExport(clus1,c("aldelval","imbensbinary","R2yzero","R2wzero","stvalnt","xmatt","W","Y"))
  clusterCall(clus1,source,file="sensitivity2MLE.R")
  
  parsapply <- function(cluster,...) {
    parSapply(cluster,...)
  }
} else {
  parsapply <- function(cluster,...) {
    sapply(...)
  }
}

system.time(
            rsc <- parsapply(clus1,1:ncomb,rsquare)
            )

dim(aldelval)
dim(t(rsc))

rsc.all <- cbind(aldelval,t(rsc))
colnames(rsc.all) <- c("alpha","delta","Tau","R2p","R2wp","R2y")
save (rsc.all,file="alpha1rsc.RData") ## saving it just in case
#Tau <- round(rsc.all[,"Tau"],1) #rare is the case we get exactly the Taus we want, rounding gives a good approximation

########################################################################
## Getting the values of Tau we are interested in (usually 2 sd + and -)
########################################################################


#Tau.in <- (Tau==0)
#
#rsc.in <- rsc.all[Tau.in,]
#
#
#
#
##########################################################################
### Calculating the Partial R2s for the observed covariates for comparison
##########################################################################
#
#totalR2w <- R2wzero
#totalR2y <- R2yzero
#
#partials <- function(xmatp){
#  stvalnp <- starting(xmatp) #starting values for MLE
#  imb <- imbensbinary(0,0,xmatp,stvalnp)
#  gammas <- as.matrix(imb@Gamma) #getting the gammas from MLE
#  betas <- as.matrix(c(imb@Beta,imb@Tau))
#  vxg <- var(xmatp[,-1]%*%gammas)
#  vxb <- var(cbind(xmatp[,-1],W)%*%betas)
#  exclR2wp <- vxg/(vxg + pi^2/3) #R2w (treatment) excluding cov of interest
#                                        #exclR2wp
#  parR2w <- (totalR2w - exclR2wp)/(1 - totalR2w) #partial R2
#
#  exclR2yp <- vxb/(vxb+pi^2/3)
#  parR2 <- (totalR2y - exclR2yp)/(1 - totalR2y) #R2 ols(outcome)excluding cov of interest
#
#  par <- cbind(parR2w,parR2)
#  par
#}
#
#X1R2p <- partials(xmatt[,-2])#xmat excluding X1
#X2R2p <- partials(xmatt[,-3])#xmat excluding X2
#X3R2p <- partials(xmatt[,-4]) #xmat excuding X3
#             
#
#covp <- rbind(X1R2p,X2R2p,X3R2p)
#
#
#
#
#
##########
### Graphs
##########
#
#ob1<-rsc.in[,4]
#ob2<-rsc.in[,5]
#
#save(ob1, file="beta8R2p.RData")
#
#save(ob2, file="beta8R2wp.RData")
#
#pdf(file="beta8graph.pdf")
#plot(rsc.in[,5],rsc.in[,4],xlim=c(0,1),xlab=expression("Partial"~ R^2 ~ "Assignment"),ylim=c(0,1),ylab=expression("Partial" ~  R^2 ~ "Outcome"))
#points(covp[,1],covp[,2],pch="+")
#graphics.off()

