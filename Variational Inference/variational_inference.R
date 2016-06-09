# if(!exists("AdmixProportion", mode="function")) source("AdMixProportion.R")
# if(!exists("AlleleFrequency", mode="function")) source("AlleleFrequency.R")
# if(!exists("marginal.likelihood", mode="function")) source("MarginalLikelihood.R")

source("AdMixProportion.R")
source("AlleleFrequency.R")
source("MarginalLikelihood.R")

infer_population_structure <- function(G, K, thres){
  total.time = proc.time()
  N = dim(G)[1]
  L = dim(G)[2]
  itertime = proc.time()
  Estart = -Inf
  
  for(restart in 1:5){
    ap = AdmixProportion(N,K)
    af = AlleleFrequency(L,K)
    
    ap = updateAPParam(ap, G, af)
    af = updateAFParam(af, G, ap)
    
    E = marginal.likelihood(G, ap, af)
    
    cat(sprintf("Marginal likelihood with initialization (%d) = %.10f\n",restart,E))
    
    if(E > Estart){
      Estart = E
      apstart = ap
      afstart = af
    }
  }
  itertime = proc.time()-itertime
  E = Estart
  ap = apstart
  af = afstart
  iter = 1
  cat(sprintf("Initial Marginal_Likelihood: %.10f Iteration_Time (secs) %.3f\n",E,itertime['elapsed']))
  

  reltol = Inf
  cat(sprintf("---------------------------------------------------------------------------\n"))
  cat(sprintf("Iteration Marginal_Likelihood delta_Marginal_Likelihood Iteration_Time (secs)\n"))
  cat(sprintf("---------------------------------------------------------------------------\n"))
  
  while(abs(reltol) > thres && iter <=150){
    ap = squareUpdateAPParam(ap,G,af)
    ap = updateAPParam(ap,G,af)
    
    af = squareUpdateAFParam(af,G,ap)
    af = updateAFParam(af,G,ap)
    
    if (iter%%10==0){
      E_new = marginal.likelihood(G, ap, af)
      reltol = E_new - E
      E = E_new
      itertime = proc.time()-itertime
      
      cat(sprintf("%d \t %.10f \t %.10f \t %.3f\n",iter,E,reltol,itertime['elapsed']))
      
      itertime = proc.time()
      #pi.update_hyperparam(False)
      
    }
    iter = iter + 1
  }
  
  P = af@var.beta/(af@var.beta + af@var.gamma)
  Q = ap@var/rowSums(ap@var)
  totaltime = proc.time() - total.time
  
  cat(sprintf("---------------------------------------------------------------------------\n"))
  cat(sprintf("Marginal Likelihood = %.10f\n",E))
  cat(sprintf("Total time = %.4f seconds\n",totaltime['elapsed']))
  cat(sprintf("Total iterations = %d \n",iter))
  
  params = list("Admixture" = Q, "AlleleFreq" = P) 
  return(params)
}