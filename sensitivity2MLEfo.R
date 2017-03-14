setClass('imbensbinary',
representation(Beta='numeric', Tau='numeric', 
                Gamma='numeric', VC='matrix', VCW='matrix',Cov='numeric'))



imbensbinary <- function(alpha,delta,xmat,stvaln) {

    logl <- function(b,alpha,delta,X,W,Y) {
        tau <- b[length(b)]
        gamma <- b[1:ncol(X)]
        beta <- b[(ncol(X)+1):(length(b)-1)]
        llik <- log(
            0.5*(((exp(X%*%beta+tau*W))^Y)/(1+exp(X%*%beta+tau*W))
               *
               ((exp(X%*%gamma))^W)/(1+exp(X%*%gamma)))
            +
            0.5*(((exp(X%*%beta+tau*W+delta))^Y)/(1+exp(X%*%beta+tau*W+delta))
               *
               ((exp(X%*%gamma+alpha))^W)/(1+exp(X%*%gamma+alpha)))
            )
    sum(llik)
    }

    imbensbinary.mle <- optim(stvaln,logl,hessian=F,method="BFGS",
    control=list(fnscale=-1,trace=1,maxit=2500,reltol=1e-17),
    alpha=alpha, delta=delta,X=xmat,W=W,Y=Y)
    
    tau <- imbensbinary.mle$par[length(stvaln)]
    beta <- imbensbinary.mle$par[(ncol(xmat)+2):(length(stvaln)-1)]
    gamma <- imbensbinary.mle$par[2:ncol(xmat)]
    gammac <- imbensbinary.mle$par[1]
    vc <- as.matrix(var(xmat[,-1]))
    vw <- cbind(xmat[,-1],W)
    vcw <- as.matrix(var(vw))

    stvaln <- imbensbinary.mle$par
    
    X<-as.matrix(xmat[,2:ncol(xmat)])
    
    #I've changed these coefficients!!!  I've added gamma here!!!
    #mod <- glm(W~X,family=binomial(link=logit))
    coef <- as.matrix(gamma)
    cont <- gammac[1]
    EUW <- .5*mean((exp(cont+X%*%coef+alpha)/(1+exp(cont+X%*%coef+alpha))))
    estcov<-EUW-.5*mean(W)  
      
    result <- new('imbensbinary',Tau=tau,Beta=beta,Gamma=gamma,VC=vc,VCW=vcw,Cov=estcov)
    class(result) <- 'imbensbinary'
    result
}

setMethod('summary', signature(object='imbensbinary'),
    definition=function(object, ...){
    Tau <- object@Tau
    table <- cbind(Tau)
    colnames(table) <- c('Tau')
    print(table)
    }
)