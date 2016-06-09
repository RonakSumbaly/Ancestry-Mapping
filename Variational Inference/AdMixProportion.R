setClass(
  "AdmixProportion",
  
  # Define the slots
  representation (
    N = "numeric",
    K = "numeric",
    alpha = "matrix",
    var = "matrix",
    xi = "matrix",
    oldvar = "list"
  )
)

AdmixProportion <- function(N, K){
  alpha = 1/K * matrix(c(1),nrow=1,ncol=K)
  var = matrix(c(1),nrow=N,ncol=K) + 0.1*matrix(runif(N*K), ncol=K)
  xi = exp(digamma(var) - digamma(rowSums(var)))
  oldvar = list()
  new("AdmixProportion",N=N,K=K,alpha=alpha,var=var,xi=xi,oldvar=oldvar)
}

get.list.tail <- function(l){
  return(l[[length(l)]])
}

makecopy.AP <- function(obj)
          {
            newObj <- AdmixProportion(
              N= obj@N,
              K=obj@K,
              alpha = obj@alpha,
              var = obj@var,
              xi = obj@xi
            )
            return(newObj)
          }

setGeneric(name="updateAPParam",
           def=function(obj,G,af)
           {
             standardGeneric("updateAPParam")
           }
)

setMethod(f="updateAPParam",
          signature="AdmixProportion",
          definition=function(obj,G, af)
          {
            #print("inside updateAPParam")
            obj@var = matrix(c(0),nrow=obj@N,ncol=obj@K)
            
            for (n in 1:obj@N) {
              # genotype = G[n,]
              # 
              # normbeta = rep(0, af@L)
              # normgamma = rep(0, af@L)

              
              # normbeta = rowSums(af@zetabeta * obj@xi[n,])
              # normgamma = rowSums(af@zetagamma * obj@xi[n,])
              # 
              # sum_a = colSums( af@zetagamma*(2-genotype) / normgamma)
              # sum_b = colSums( af@zetabeta*genotype / normbeta)
              # 
              # obj@var[n,] = obj@var[n,] + (sum_a + sum_b) * obj@xi[n,]
              
              for (l in 1:af@L) {
                normbeta = 0.0
                normgamma = 0.0
                genotype = G[n,l]
                #  compute normalization
                #for (k in 1:obj@K) {
                #  normbeta = normbeta + af@zetabeta[l,k] * obj@xi[n,k]
                #  normgamma = normgamma + af@zetagamma[l,k] * obj@xi[n,k]
                #}
                normbeta = sum(af@zetabeta[l,] * obj@xi[n,])
                normgamma = sum(af@zetagamma[l,] * obj@xi[n,])
                # print(normbeta)
                # loop over populations
                #for (k in 1:obj@K) {
                #compute new estimate of variational parameters
                #   obj@var[n,k] = obj@var[n,k] + (((2-genotype) * af@zetagamma[l,k] / normgamma ) + (genotype * af@zetabeta[l,k] /  normbeta) ) * obj@xi[n,k]
                #}
                 obj@var[n,] = obj@var[n,] + (((2-genotype) * af@zetagamma[l,] / normgamma ) + (genotype * af@zetabeta[l,] /  normbeta) ) * obj@xi[n,]

               }
            }
            
            if(anyNA(obj@var)){
              obj@var = get.list.tail(obj@oldvar)
            }else{
              obj@var = rep(obj@alpha, obj@N) + obj@var
              obj@xi = exp(digamma(obj@var) - digamma(rowSums(obj@var)))
            }
            return(obj)
          }
)
setGeneric(name="squareUpdateAPParam",
           def=function(obj,G,af)
           {
             standardGeneric("squareUpdateAPParam")
           }
)

setMethod(f="squareUpdateAPParam",
          signature="AdmixProportion",
          definition=function(obj,G,af)
          {
            #print("inside squareUpdateAPParam")
            obj@oldvar = list(obj@var)
            # take two update steps
            for(step in 1:2){
              obj = updateAPParam(obj, G, af)
              obj@oldvar = c(obj@oldvar, list(obj@var))
            }
              
            R = obj@oldvar[[2]] - obj@oldvar[[1]]
            V = obj@oldvar[[3]] - obj@oldvar[[2]] - R
            a = -1 * sqrt(sum(R*R)/sum(V*V))
            
            if(a>-1){
              a = -1
            }
              
            # given two update steps, compute an optimal step that achieves
            # a better marginal likelihood than the best of the two steps.
            a_ok = FALSE
            while(!a_ok){
              obj@var = (1+a)^2*obj@oldvar[[1]] - 2*a*(1+a)*obj@oldvar[[2]] + a^2*obj@oldvar[[3]]
              if (all(obj@var > 0)){
                a_ok = TRUE
              }else{
                a = (a-1)/2.
                if(abs(a+1)< exp(-4)){
                  a = -1
                 }
              }
            }
              
            # if obj accelerated step fails for some reason, stick with the first non-accelerated step.
            if(anyNA(obj@var)){
              obj@var = obj@oldvar[[2]]
            }
            obj@xi = exp(digamma(obj@var) - digamma(rowSums(obj@var)))
            return(obj)
          }
)