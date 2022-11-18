# Statistical Programming group project 4 code 
# Group 35: 
# Henry Blackwell (s2451994) comments and pusedocode/logic, Peijie Zeng (s2332799) , Jing Pan (s2312688)
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
## Overview:
## Function 'newt' is a function which implements Newton's method to a given function
## with its gradient, Hessian, function and derivatives. 'newt' also
## takes some inputs includes initial values, convergence tolerance, maximum
## iteration time and other coefficients to run the optimisation process.
##
## Inputs:
## theta: Initial points to start the process, could be a list of values
## func: a function with parameters of the objective function
## grad: a function with parameters of the gradient of the objective function
## hess: a function of Hessian matrix of the objective function, if not given,
## newt will get the approximation with finite difference of the gradient

## ... , arguments passed to func, grad and hess
## tol: the convergence tolerance used to check whether find the optimum
## fscale: a rough estimate of the magnitude of 'func' near the optimum
## mixit: the maximum number of Newton iterations to try
## max.half: the maximum number of times a step should be halved before
## concluding that the step has failed to improve the objective.
## eps: the finite difference intervals used to approximate the Hessian matrix
##
## Outpus:
## f: the value of the objective function at the minimum
## theta: value of the parameters at the minimum
## iter: the number of iterations taken to reach the minimum
## g: the gradient vector at the minimum
## Hi: the inverse of the Hessian matrix at the minimum



{
  dim=length(theta) # gets the number of parameters
  iter=1 # intiliases the interations to 1
  # while interations are less then the user defined maximum, perform newton opti
  while (iter<=maxit) 
  {
    func_val0=func(theta,...) # set func_val0 equal to input
    grad_val=grad(theta,...) # grad_val equal to input
    # checks if the intial value of the function isn't finite
    if (is.finite(func_val0)!=TRUE) 
    {
      stop('Initial point doesn\'t fit') # stop and print an error message 
    }
    if (is.null(hess)) # if we aren't given a hessian matrix, use finite differencing
    {
      hess_val=matrix(0,dim,dim) 
      for (i in 1:dim) # these loops are to go over and "calulate the partials for each 
                       # parameter in theta, i loops over row and j over columns
                       #  only calulate the upper triangler, as it's more effecient
      {
        for (j in i:dim)
        {
          ei=rep(0,dim) # finite increasement in the ith dimension, delta*x*i 
          ei[i]=1
          ej=rep(0,dim) # ''' x*j 
          ej[j]=1
          grad_pi=grad(theta+ei*eps,...)  # fin the finite diff intergrals 
          grad_mi=grad(theta-ei*eps,...)   
          grad_pj=grad(theta+ej*eps,...)
          grad_mj=grad(theta-ej*eps,...)
          # the main equation for calulating the hessian at each position 
          hess_val[i,j]=(grad_pi[j]-grad_mi[j]+grad_pj[i]-grad_mj[i])/(4*eps)
        }
      }
      hess_val=(hess_val+t(hess_val))/2 # makes sure it's symmetric   
    }
    else # if we are given the hessian then just use that obvoiusly 
    {
      hess_val=hess(theta,...)  
    }
    # take the egin values of the hessain
    diag_elements=eigen(hess_val)$values 
    # optimum checking (if we really are at a minimum at not saddle point etc)
    # checks each value against the given tolerance 
    if (all(abs(grad_val)<tol*(abs(func_val0)+fscale)))
    {
      if (any(diag_elements<=0)) #  the hessain ins't postive definite
      {
        stop('failed to converge') # print error message as it can't be a minimum
        # maybe just a saddle etc etc
        break
      }
      else # we are at a minimum 
      {
        break # hence break becuasse we have completed our objective
      }
    }
    iter=iter+1 
    if (iter==(maxit+1)) # if we can't find the optimal value in maxit trys (+1)
    {
      stop('can\'t find the optimum in given steps') # then display an error 
    }
    # if the given hessian isn't poistive definite (pd), perturb until it is
    perturb_m=0 # this is a constant which actullay used to increase the perturbment multiplier 
    while (any(diag_elements<=0)) # while egin values aren't pd
    {
      # adding an idenitity matrix to achive pd
      hess_val=hess_val+10^(perturb_m-6)*diag(1,dim) 
      # extract the egin values from the new perturbed matix
      diag_elements=eigen(hess_val)$values 
      perturb_m=perturb_m+1 # increment this perturbment multiplier 
    }
    # creating a step length using cholsky decom
    delta=-chol2inv(chol(hess_val))%*%grad_val 
    
    # check whether to halve the step
    for (i in 1:(max.half+1)) 
    {
      # delta equals to the of negative of ( gradient * inverse hessian )
      # delata is the intiial choose of the step length 
      theta_test=theta+delta
      fun_val1=func(theta_test,...) # update the function value 
      # checks if we can ommite the max half limit, if true then return error
      # we also checks if the function is finite
      if (i==max.half+1 || is.finite(fun_val1)!=TRUE) 
      {
        stop('step failed to optimise the function') # display the given probelm
      }
      if (fun_val1>=func_val0) # if the step doesn't optimises then half the step
      {
        delta=delta/2 # half the step
      }
      else
      {
        # if we do optimise the function then, set theta as the optimal parameter value 
        theta=theta_test 
        break
      }
    }
  }
  # return, the value of the function at the optimal point, optimal parameter values
  # and the gradient, inverse of the hessian at the optimal point and the number of iterations
  result=list(f=func_val0,
              theta=theta,iter=iter,
              g=grad_val,
              Hi=chol2inv(chol(hess_val)))
  return(result)
}