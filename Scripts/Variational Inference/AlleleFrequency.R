setClass(
  "AlleleFrequency",
  
  # Define the slots
  representation (
    L = "numeric",
    K = "numeric",
    beta = "matrix",
    gamma = "matrix",
    var.beta = "matrix",
    var.gamma = "matrix", 
    zetabeta = "matrix",
    zetagamma = "matrix",
    oldvar.beta = "list",
    oldvar.gamma = "list"
  )
)

AlleleFrequency <- function(L, K){
  beta = matrix(c(1),nrow=L,ncol=K)
  gamma = matrix(c(1),nrow=L,ncol=K)
  var.beta = matrix(c(1),nrow=L,ncol=K) + 0.1*matrix(runif(L*K), ncol=K)
  var.gamma = matrix(c(10),nrow=L,ncol=K) + 0.1*matrix(runif(L*K), ncol=K)
  zetabeta = exp(digamma(var.beta) - digamma(var.beta + var.gamma))
  zetagamma = exp(digamma(var.gamma) - digamma(var.beta + var.gamma))
  oldvar.beta = list()
  oldvar.gamma = list()
  new("AlleleFrequency",L=L,K=K,beta=beta, gamma = gamma, var.beta = var.beta, var.gamma = var.gamma, 
      zetabeta = zetabeta, zetagamma= zetagamma, oldvar.beta=oldvar.beta, oldvar.gamma=oldvar.gamma)
}

makecopy.AF <- function(obj)
{
  newObj <- AlleleFrequency(
    L = obj@L,
    K = obj@K,
    beta = obj@beta,
    gamma = obj@gamma,
    var.beta = obj@var.beta,
    var.gamma = obj@var.gamma,
    zetabeta = obj@zetabeta,
    zetagamma = obj@zetagamma,
    oldvar.beta = obj@oldvar.beta,
    oldvar.gamma = obj@oldvar.gamma
  )
  return(newObj)
}

setGeneric(name="updateAFParam",
           def=function(obj,G,ap)
           {
             standardGeneric("updateAFParam")
           }
)

setMethod(f="updateAFParam",
          signature="AlleleFrequency",
          definition=function(obj,G, ap)
          {
            #print("inside updateAFParam")
            obj@var.beta = matrix(c(0),nrow=obj@L,ncol=obj@K)
            obj@var.gamma = matrix(c(0),nrow=obj@L,ncol=obj@K)
            
            var.beta.tmp <- rep(0, obj@K)
            var.gamma.tmp <- rep(0, obj@K)
            
            # loop over loci
            for (l in 1:obj@L) {
              
              #for (k in 1:obj@K) {
              #  
              #  var.beta.tmp[k] = 0.0
              #  var.gamma.tmp[k] = 0.0
              #}
              
              var.beta.tmp <- rep(0, obj@K)
              var.gamma.tmp <- rep(0, obj@K)
              
              genotype = G[,l]
              
              theta.beta.sum = rep(0, ap@N)
              theta.gamma.sum = rep(0, ap@N)
              
              # theta.beta.sum = rowSums(t(t(ap@xi) * obj@zetabeta[l,]))
              # theta.gamma.sum = rowSums(t(t(ap@xi) * obj@zetagamma[l,]))
              
              theta.beta.sum = rowSums(ap@xi * obj@zetabeta[l,])
              theta.gamma.sum = rowSums(ap@xi * obj@zetagamma[l,])
              
              # var.beta.tmp = colSums( t(t(t(t(ap@xi)*genotype)) / theta.beta.sum))
              # var.gamma.tmp = colSums( t(t(t(t(ap@xi)*(2-genotype))) / theta.gamma.sum))
              
              var.beta.tmp = colSums( ap@xi*genotype / theta.beta.sum)
              var.gamma.tmp = colSums( ap@xi*(2-genotype) / theta.gamma.sum)
              
              # loop over samples
              #for (n in 1:ap@N) {
                
                # genotype = G[n,l]
                
                
                # compute xi*zeta_{beta,gamma}
                # theta.beta.sum = 0.0
                # theta.gamma.sum = 0.0
                # 
                #for (k in 1:obj@K) {
                #  theta.beta.sum = theta.beta.sum + ap@xi[n,k] * obj@zetabeta[l,k]
                #  theta.gamma.sum = theta.gamma.sum + ap@xi[n,k] * obj@zetagamma[l,k]
                #}
                #theta.beta.sum = sum(ap@xi[n,] * obj@zetabeta[l,])
                #theta.gamma.sum = sum(ap@xi[n,] * obj@zetagamma[l,])
                
                # increment var_{beta,gamma}_tmp
                #for (k in 1:obj@K) {
                #  var.beta.tmp[k] = var.beta.tmp[k] +  genotype * ap@xi[n,k] / theta.beta.sum
                #  var.gamma.tmp[k] = var.gamma.tmp[k] + (2-genotype) * ap@xi[n,k] / theta.gamma.sum
                #}
                # var.beta.tmp = var.beta.tmp +  genotype * ap@xi[n,] / theta.beta.sum
                # var.gamma.tmp = var.beta.tmp + (2-genotype) * ap@xi[n,] / theta.gamma.sum
              #}
              
              # compute var_{beta,gamma}
              #for (k in 1:obj@K) {
              #  obj@var.beta[l,k] = obj@beta[l,k] + obj@zetabeta[l,k] * var.beta.tmp[k]
              #  obj@var.gamma[l,k] = obj@gamma[l,k] + obj@zetagamma[l,k] * var.gamma.tmp[k]
              #}
              
              obj@var.beta[l,] = obj@beta[l,] + obj@zetabeta[l,] * var.beta.tmp
              obj@var.gamma[l,] = obj@gamma[l,] + obj@zetagamma[l,] * var.gamma.tmp
            }
            if(anyNA(obj@var.beta)){
              obj@var.beta = get.list.tail(obj@oldvar.beta)
            }
            if(anyNA(obj@var.gamma)){
              obj@var.gamma = get.list.tail(obj@oldvar.gamma)
            }
            
            obj@zetabeta = exp(digamma(obj@var.beta) - digamma(obj@var.beta + obj@var.gamma))
            obj@zetagamma = exp(digamma(obj@var.gamma) - digamma(obj@var.beta + obj@var.gamma))
            
            return(obj)
          }
)
setGeneric(name="squareUpdateAFParam",
           def=function(obj,G,ap)
           {
             standardGeneric("squareUpdateAFParam")
           }
)

setMethod(f="squareUpdateAFParam",
          signature="AlleleFrequency",
          definition=function(obj,G,ap)
          {
            #print("inside squareUpdateAFParam")
            obj@oldvar.beta = list(obj@var.beta)
            obj@oldvar.gamma = list(obj@var.gamma)
            
            # take two update steps
            for(step in 1:2){
              obj = updateAFParam(obj,G,ap)
              obj@oldvar.beta = c(obj@oldvar.beta, list(obj@var.beta))
              obj@oldvar.gamma = c(obj@oldvar.gamma, list(obj@var.gamma))
            }
              
            R.beta = obj@oldvar.beta[[2]] - obj@oldvar.beta[[1]]
            R.gamma = obj@oldvar.gamma[[2]] - obj@oldvar.gamma[[1]]
            V.beta = obj@oldvar.beta[[3]] - obj@oldvar.beta[[2]] - R.beta
            V.gamma = obj@oldvar.gamma[[3]] - obj@oldvar.gamma[[2]] - R.gamma
            
            a = -1* sqrt((sum(R.beta*R.beta) + sum(R.gamma*R.gamma)) / (sum(V.beta*V.beta) + sum(V.gamma*V.gamma)))
            
            if(a>-1){
              a = -1 
            }
              
            
            # given two update steps, compute an optimal step that achieves
            # a better marginal likelihood than the best of the two steps.
            a_ok = FALSE
            while(!a_ok){
              obj@var.beta = (1+a)^2*obj@oldvar.beta[[1]] - 2*a*(1+a)*obj@oldvar.beta[[2]] + a^2*obj@oldvar.beta[[3]]
              obj@var.gamma = (1+a)^2*obj@oldvar.gamma[[1]] - 2*a*(1+a)*obj@oldvar.gamma[[2]] + a^2*obj@oldvar.gamma[[3]]
              if (all(obj@var.beta > 0) && all(obj@var.gamma > 0)){
                a_ok = TRUE
              }else{
                a = (a-1)/2
                if(abs(a+1)< exp(-4)){
                  a = -1
                }
              }
            }
            
            # if this accelerated step fails for some reason, stick with the first non-accelerated step.
            if(anyNA(obj@var.beta) || anyNA(obj@var.gamma)){
              obj@var.beta = obj@oldvar.beta[[2]]
              obj@var.gamma = obj@oldvar.gamma[[2]]
            }
              
            obj@zetabeta = exp(digamma(obj@var.beta) - digamma(obj@var.beta + obj@var.gamma))
            obj@zetagamma = exp(digamma(obj@var.gamma) - digamma(obj@var.beta + obj@var.gamma))
            
            return(obj)
          }
)
