# Statistical Programming group project 4 code 
# Group 35: 
# Henry Blackwell (s2451994), Peijie Zeng (s2332799) , Jing Pan (s2312688)
# https://github.com/zengens/group35_proj4
# Peijie Zeng: 
#
#
# READ ME: 
# For function description, outline and overview (using one: #)
# line-by-line comments will above the corresponding code (using two: ## )
# Multiple lines of comments may be use for explanation
#

newt=function(theta,func,grad,hess,...,tol,fscale,maxit,max.half,eps)
{
  dim=length(theta)
  iter=1
  while (iter<=maxit)
  {
    func_val0=func(theta,...)
    grad_val=grad(theta,...)
    if (is.finite(func_val0)!=TRUE)
    {
      stop('Initial point doesn\'t fit')
    }
    if (is.null(hess))
    {
      hess_val=matrix(0,dim,dim)
      for (i in 1:dim)
      {
        for (j in i:dim)
        {
          ei=rep(0,dim)
          ei[i]=1
          ej=rep(0,dim)
          ej[j]=1
          grad_pi=grad(theta+ei*eps,...)
          grad_mi=grad(theta-ei*eps,...)
          grad_pj=grad(theta+ej*eps,...)
          grad_mj=grad(theta-ej*eps,...)
          hess_val[i,j]=(grad_pi[j]-grad_mi[j]+grad_pj[i]-grad_mj[i])/(4*eps)
        }
      }
      hess_val=(hess_val+t(hess_val))/2
    }
    else 
    {
      hess_val=hess(theta,...)
    }
    diag_elements=eigen(hess_val)$values
    # optimum checking
    if (all(abs(grad_val)<tol*(abs(func_val0)+fscale)))
    {
      if (any(diag_elements<=0))
      {
        stop('failed to converge')
        break
      }
      else
      {
        break
      }
    }
    iter=iter+1
    if (iter==(maxit+1))
    {
      stop('can\'t find the optimum in given steps')
    }
    perturb_m=0
    while (any(diag_elements<=0))
    {
      perturb_m=perturb_m+1
      hess_val=hess_val+10^(perturb_m-6)*diag(1,dim)
      diag_elements=eigen(hess_val)$values
    }
    delta=-chol2inv(chol(hess_val))%*%grad_val
    # check whether to halve the step
    for (i in 1:(max.half+1))
    {
      
      theta_test=theta+delta
      fun_val1=func(theta_test,...)
      if (i==max.half+1 || is.finite(fun_val1)!=TRUE)
      {
        stop('step failed to optimise the function')
      }
      if (fun_val1>=func_val0)
      {
        delta=delta/2
      }
      else
      {
        theta=theta_test
        break
      }
    }
  }
  result=list(f=func_val0,
              theta=theta,iter=iter,
              g=grad_val,
              Hi=chol2inv(chol(hess_val)))
  return(result)
}