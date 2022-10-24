


###############################################################################################
###############################################################################################
# 
# Load functions
#
###############################################################################################
###############################################################################################


#----------------------------------------------------------
# 
# Structure of lambda: 6 rows, m-1 columns
#
#----------------------------------------------------------
#
# lambda = [
# 	alpha	. . . alpha;
#	p1	. . . p1;
#	q1	. . . q1;
#	p2	. . . pm;
#	q2	. . . qm;
#	w	. . . wm]
#
#----------------------------------------------------------



###############################################################################################
# Auxiliary functions
###############################################################################################


#-------------------------------------------------------------
#symbolic Higher derivatives:
#-------------------------------------------------------------

DD = function(expr, name, order = 1) {
  if(order < 1){ stop("'order' must be >= 1") }
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}


#-------------------------------------------------------------
# Higher derivatives of u  at zero:
#-------------------------------------------------------------

der_u0_SingP = function(lambda, n){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r= (-q2*w)/q1
  
  u = expression((p1*exp(sicma*x) + q1*exp((sicma-eta1)*x))^r)
  x=0; 
  u0 = eval(u)
  
  for (i in 1:(n-1)) {
    g=DD(u,"x",i)
    x=0
    u0=c(u0,eval(g))
  } 
  return(u0) 
}


#-------------------------------------------------------------
# Build matrix hat_B
#-------------------------------------------------------------

Build_B_SingP = function(lambda, n){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r= (-q2*w)/q1
  
  ## Matrix B
  
  hat_B = matrix(0,n,n)
  i = 0:(n-1)
  
  for(j in 1:n){
    hat_B[,j] = (j*eta1 + (sicma - eta1)*r)^i
  }
  
  
  return(hat_B)
  
}


#-------------------------------------------------------------
# Hager inverse matrix
#-------------------------------------------------------------

H_inv_B_SingP = function(lambda, n){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
  b0 = rep(bj0,n)^(0:(n-1)) 
  
  ## Matrix hat_B
  hat_B = Build_B_SingP(lambda,n)
  
  
  if(n==2){
    
    inv_hat_B = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
    
  }else{
    
    ## Partitions of matrix hat_B 
    #hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
    
    ## Partitions of the inverse of matrix hat_B  
    A = hat_B3 - hat_B2%*%(inv_hat_B0%*%hat_B1)
    tilde_invB0 = inv_hat_B0 + (inv_hat_B0%*%hat_B1)%*%solve(A)%*%(hat_B2%*%inv_hat_B0)
    tilde_invB1 = inv_hat_B0%*%hat_B1%*%solve(A)
    tilde_invB2 = solve(A)%*%hat_B2%*%inv_hat_B0
    
    ## Inverse
    M1 = rbind(tilde_invB0,- tilde_invB2)
    M2 = rbind(-tilde_invB1, solve(A))
    inv_hat_B = cbind(M1,M2)
  }
  return(inv_hat_B)
}

###############################################################################################
# S and Delta functions from Lemma 6
###############################################################################################

#-------------------------------------------------------------
# Function Delta1 (from Lemma 6)
#-------------------------------------------------------------

Delta1 = function(lambda,n,t){
  
  if(n<=2){
    
    Del1 = rep(0,length(t))
    
  }else{
    
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    
    bj0 = (sicma - eta1)*r
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A = hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0 = (solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1 = solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2 = solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta1
    M1 = rbind(tilde_invB0, -tilde_invB2)
    M2 = rbind(-tilde_invB1, solve(A))
    M = cbind(M1,M2)
    
    Del1=NULL
    for (k in 1:length(t)) {
      #exp_bt = exp(b_j*t[k])
      Del1 = c(Del1,as.numeric(t(exp(b_j*t[k]))%*%(M%*%b0)))
    }
    
    
  }
  
  return(Del1)
}

#-------------------------------------------------------------
# S_1 (from Lemma 6)
#-------------------------------------------------------------

S_1 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  a=NULL
  for (k in 1:length(t)) {
	  
    #exp_bt=exp(b_j*t[k])
    a=c(a,as.numeric( t(exp(b_j*t[k]))%*%(inv_hat_B0%*%b0)))
    
  }
  
  
  
  return(a)
}


#-------------------------------------------------------------
# Function Delta2 (from Lemma 6)
#-------------------------------------------------------------

Delta2 = function(lambda,n,t){
  
  if(n<=2){
    Del2= rep(0,length(t))
  }else{
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    t1 = log(q1/p1)/eta1
    bj0 = (sicma - eta1)*r
    
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A = hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0 = (solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1 = solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2 = solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta1
    M1 = rbind(tilde_invB0,- tilde_invB2)
    M2 = rbind(-tilde_invB1, solve(A))
    M = cbind(M1,M2)
    
    Del2=NULL
    for (k in 1:length(t)) {
      Del2=c(Del2,as.numeric(t(exp(b_j*t[k])/b_j)%*%(M%*%b0)))
    }
  }
  
  
  return(Del2)
}

#-------------------------------------------------------------
# S_2 (from Lemma 6)
#-------------------------------------------------------------

S_2 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  a = NULL
  for (k in 1:length(t)) {
    exp_bt = exp(b_j*t[k])/b_j
    a = c(a,as.numeric(t(exp(b_j*t[k])/b_j)%*%(inv_hat_B0%*%b0)))
  }
  
  # rm(p1, q1, p2, q2, w, sicma, r, eta1, eta2, bj0, b_j, b0, lambda, t1, exp_bt)
  
  return(a)
}


#-------------------------------------------------------------
# Function Delta3 (from Lemma 6)
#-------------------------------------------------------------

Delta3 = function(lambda,n,t){
  
  if(n<=2){
    
    Del3=rep(0,length(t))
    
  }else{
    
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    
    bj0 = (sicma - eta1)*r
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## values at point 0 for n derivatives of the function u 
    u0 = der_u0_SingP(lambda,n) 
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A = hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0 = (solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1 = solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2 = solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta3
    M1 =rbind(tilde_invB0,-tilde_invB2)
    M2 =rbind(-tilde_invB1,solve(A))
    M = cbind(M1,M2)
    
    Del3=NULL
    for (k in 1:length(t)) {
      exp_bt = exp(b_j*t[k])
      Del3 = c(Del3,as.numeric(t(exp_bt)%*%(M%*%u0)))
    }
    
    
  }
  
  return(Del3)
}

#-------------------------------------------------------------
# S_3 (from Lemma 6)
#-------------------------------------------------------------

S_3 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  ## values at point 0 for n derivatives of the function u 
  u0 = der_u0_SingP(lambda,2)
  
  a = NULL
  for (k in 1:length(t)) {
    a = c(a,as.numeric(t(exp(b_j*t[k]))%*%(inv_hat_B0%*%u0) ))
  }
  
  
  return(a)
  
}


#-------------------------------------------------------------
# Function Delta4 (from Lemma 6)
#-------------------------------------------------------------

Delta4 = function(lambda,n,t){
  
  if(n<=2){
    Del4= rep(0,length(t))
  }else{
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    t1 = log(q1/p1)/eta1
    bj0 = (sicma - eta1)*r
    
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## values at point 0 for n derivatives of the function u 
    u0 = der_u0_SingP(lambda,n)
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A= hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0=(solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1=solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2=solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta3
    M1 = rbind(tilde_invB0,-tilde_invB2)
    M2 = rbind(-tilde_invB1,solve(A))
    M = cbind(M1,M2)
    
    Del4=NULL
    for (k in 1:length(t)) {
      Del4=c(Del4,as.numeric(t(exp(b_j*t[k])/b_j)%*%(M%*%u0)))
    }
  }
  
  
  return(Del4)
}

#-------------------------------------------------------------
# S_4 (from Lemma 6)
#-------------------------------------------------------------

S_4 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  ## values at point 0 for n derivatives of the function u 
  u0 = der_u0_SingP(lambda,2)
  
  a = NULL
  for (k in 1:length(t)) {
    a = c(a,as.numeric(t(exp(b_j*t[k])/b_j)%*%(inv_hat_B0%*%u0)))
  }
  
  return(a)
  
}


#-------------------------------------------------------------
# Function Delta5 (from Lemma 6)
#-------------------------------------------------------------

Delta5 = function(lambda,n,t){
  
  if(n<=2){
    Del5= rep(0,length(t))
  }else{
    
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    t1 = log(q1/p1)/eta1
    bj0 = (sicma - eta1)*r
    
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A= hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0=(solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1=solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2=solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta5
    M1 = rbind(tilde_invB0,-tilde_invB2)
    M2 = rbind(-tilde_invB1,solve(A))
    M = cbind(M1,M2)
    
    Del5=NULL
    for (k in 1:length(t)) {
      Del5=c(Del5,as.numeric(t(b_j*exp(b_j*t[k]))%*%(M%*%b0)))
    }
  }
  
  
  return(Del5)
}

#-------------------------------------------------------------
# S_5 (from Lemma 6)
#-------------------------------------------------------------

S_5 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  a = NULL
  for (k in 1:length(t)) {
    a = c(a,as.numeric(t(b_j*exp(b_j*t[k]))%*%(inv_hat_B0%*%b0)))
  }
  
  return(a)
  
}


#-------------------------------------------------------------
# Function Delta6 (from Lemma 6)
#-------------------------------------------------------------

Delta6 = function(lambda,n,t){
  
	if(n<=2){
		Del6= rep(0,length(t))
	}else{
    
    p1 = lambda[1]
    q1 = lambda[2]
    p2 = lambda[3]
    q2 = lambda[4]
    w = lambda[5]
    
    w = max(1.0e-15, w)
    
    eta2 = p2+q2
    eta1 = p1+q1
    
    sicma = (eta2*q1)/(q2*w)
    r = (-q2*w)/q1
    
    bj0 = (sicma - eta1)*r
    b_j = rep(eta1,n)*(1:n) + rep(bj0,n)
    b0 = rep(bj0,n)^(0:(n-1)) 
    
    ## Matrix hat_B
    hat_B = Build_B_SingP(lambda,n)
    
    ## values at point 0 for n derivatives of the function u 
    u0 = der_u0_SingP(lambda,n) 
    
    ## Partitions of matrix hat_B 
    hat_B0 = hat_B[1:2,1:2]
    hat_B1 = hat_B[1:2,3:n]
    hat_B2 = hat_B[3:n,1:2]
    hat_B3 = hat_B[3:n,3:n]
    
    ## Partitions of the inverse of matrix hat_B  
    A= hat_B3 - hat_B2%*%(solve(hat_B0)%*%hat_B1)
    tilde_invB0=(solve(hat_B0)%*%hat_B1)%*%solve(A)%*%(hat_B2%*%solve(hat_B0))
    tilde_invB1=solve(hat_B0)%*%hat_B1%*%solve(A)
    tilde_invB2=solve(A)%*%hat_B2%*%solve(hat_B0)
    
    ## value of delta3
    M1=rbind(tilde_invB0,-tilde_invB2)
    M2=rbind(-tilde_invB1,solve(A))
    M= cbind(M1,M2)
    
    Del6=NULL
    for (k in 1:length(t)) {
      Del6=c(Del6,as.numeric(t(b_j*exp(b_j*t[k]))%*%(M%*%u0)))
    }
  }
  
  
  return(Del6)
  
}

#-------------------------------------------------------------
# S_6 (from Lemma 6)
#-------------------------------------------------------------

S_6 = function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  #w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  
  bj0 = (sicma - eta1)*r
  b_j = rep(eta1,2)*(1:2) + rep(bj0,2)
  b0 = rep(bj0,2)^(0:1) 
  
  ## inverse matrix hat_B0 
  inv_hat_B0 = rbind(c(b_j[2],-1),c(-b_j[1],1))/eta1
  
  ## values at point 0 for n derivatives of the function u 
  u0 = der_u0_SingP(lambda,2)
  
  a = NULL
  for (k in 1:length(t)) {
    a = c(a,as.numeric(t(b_j*exp(b_j*t[k]))%*%(inv_hat_B0%*%u0)))
  }
  
  
  return(a)
}




###############################################################################################
# F function
###############################################################################################


#-------------------------------------------------------------
# Numerator in the F function
#-------------------------------------------------------------

Num_True_SingP =function(lambda,t){
  
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r= (-q2*w)/q1
  h = (p1*exp(sicma*t) + q1*exp((sicma-eta1)*t))^r
  
  
  return(h)
  
}

#-------------------------------------------------------------
# Denominator in the F function
#-------------------------------------------------------------

Den_True_SingP = function(lambda,t){
  
  #lambda = [alpha,p1,q1,p2,q2,w]
  
  p1 = lambda[1]
  q1 = lambda[2]
  eta1 = p1+q1
  p2 = lambda[3]
  q2 = lambda[4]
  eta2 = p2+q2
  w = lambda[5]
  sicma = (eta2*q1)/(q2*w)
  r= (-q2*w)/q1
  
  h = rep(0,length(t))
  for (i in 1:length(t)) {
    h[i] = q2*(1-w)*integrate(Num_True_SingP, 0, t[i], lambda = lambda)$value - Num_True_SingP(t = 0, lambda = lambda)
  }
  
  return(h)
  
  
}

#-------------------------------------------------------------
# The F function
#-------------------------------------------------------------

F_True_SingP = function(lambda,t){ 
  
  1 + Num_True_SingP(lambda,t)/Den_True_SingP(lambda,t)
  
}

#---------------------------------------------------------



#----------------------------------------------------------------------------------------
# Approximation of the numerator in the true function before t_max
# based on single point precision
#----------------------------------------------------------------------------------------

# Non normalized function

hat_u1_SingP = function(lambda,n,t,Expression = 2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  b0 = bj0^(0:(n-1))
  g_j0 = (q1/p1)^(b_j/eta1)
  
  
  inv_B = H_inv_B_SingP(lambda,n)
  u0 = der_u0_SingP(lambda,n)
  t1 = log(q1/p1)/(q1 + p1)
  
  C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (t(g_j0)%*%(inv_B%*%u0))
  C0_d = ((q1/p1)^(bj0/eta1)) - (t(g_j0)%*%(inv_B%*%b0))
  C_0  = C0_u/C0_d
  
  H = NULL
  
  for (i in 1:L_t){
    	
    if(Expression == 2){ # Delta decomposition from the expression of hat_u1 in the proof of Proposition 1
      
      C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - S_3(lambda,t1) - Delta3(lambda,n,t1)
      C0_d = ((q1/p1)^(bj0/eta1)) - S_1(lambda,t1) - Delta1(lambda,n,t1)
      C_0  = C0_u/C0_d
      
      H = c(H, as.numeric(C_0*( exp(bj0*t[i]) - S_1(lambda,t[i]) - Delta1(lambda,n,t[i]) ) + S_3(lambda,t[i]) + Delta3(lambda,n,t[i])))
      
    }
    
  }  
  
  
  return(H)
}


N1_SingP = function(lambda,n,t,Expression = 2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  
  t1 = log(q1/p1)/(q1 + p1)
  
  K_plus  = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (eta1^r)*( (2-(p1*r)/eta1) - (1-(p1*r)/eta1)*(q1/p1) )*((q1/p1)^(b_j[1]/eta1))
  K_minus = ((q1/p1)^(bj0/eta1)) - (2 - (q1/p1))*((q1/p1)^(b_j[1]/eta1))
  C0_u    = K_plus - Delta3(lambda,n,t1)
  C0_d    = K_minus- Delta1(lambda,n,t1)
  
  tilde_K_plus  = ((2*q1)^r)*((p1^2)/q1) - (eta1^r)*( p1*(2-(p1*r)/eta1) - q1*(1-(p1*r)/eta1) )
  tilde_K_minus = ((p1^2)/q1) + q1 - 2*p1
  
  H = NULL
  
  for (i in 1:L_t){
    const1 = p1*((p1/q1)^(b_j[1]/eta1))*exp(-b_j[1]*t[i])
    
    if(Expression == 2){ # Delta decomposition from the expression of hat_u1 in the proof of Proposition 1
      
      tilde_N1_plus  = exp(-eta1*t[i]) + exp(eta1*t[i]) - 2
      tilde_N1_minus = (eta1^r)*( (2-(p1*r)/eta1) - (1-(p1*r)/eta1)*exp(eta1*t[i]) )
      
      H = c(H, as.numeric( tilde_K_plus*tilde_N1_plus + tilde_K_minus*tilde_N1_minus + 
				const1*(C0_d*Delta3(lambda,n,t[i]) - C0_u* Delta1(lambda,n,t[i]) - ( exp(bj0*t[i]) - S_1(lambda,t[i]))*Delta3(lambda,n,t1) -  
				S_3(lambda,t[i])*Delta1(lambda,n,t1) ) ) )
    }

  }  
  
  return(H)
  
}



#----------------------------------------------------------------------------------------
# Approximation of the denominator in the true function before t_max
#----------------------------------------------------------------------------------------

# Non normalized function

hat_den1_SingP = function(lambda,n,t, Expression = 2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  b0 = bj0^(0:(n-1))
  g_j0 = (q1/p1)^(b_j/eta1)
  
  
  inv_B = H_inv_B_SingP(lambda,n)
  u0 = der_u0_SingP(lambda,n)
  t1 = log(q1/p1)/(q1 + p1)
  
  C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (t(g_j0)%*%(inv_B%*%u0))
  C0_d = ((q1/p1)^(bj0/eta1)) - (t(g_j0)%*%(inv_B%*%b0))
  C_0  = C0_u/C0_d
  
  H = NULL
  
  for (i in 1:L_t){
    	
    if(Expression == 2){ # Delta decomposition from the expression of hat_u1 in the proof of Proposition 1
      
      C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - S_3(lambda,t1) - Delta3(lambda,n,t1)
      C0_d = ((q1/p1)^(bj0/eta1)) - S_1(lambda,t1) - Delta1(lambda,n,t1)
      C_0  = C0_u/C0_d
      
      kappa_02 = as.numeric(C_0*( 1 - S_1(lambda,0) - Delta1(lambda,n,0) ) +  S_3(lambda,0) + Delta3(lambda,n,0)+
                              q2*(1-w)*( C_0*( (1/bj0)  - S_2(lambda,0) - Delta2(lambda,n,0) ) +  S_4(lambda,0) + Delta4(lambda,n,0) ) )   
      
      H = c(H, as.numeric( q2*(1-w)*( C_0*( (exp(bj0*t[i])/bj0)  - S_2(lambda,t[i]) - Delta2(lambda,n,t[i]) ) +  
                                        S_4(lambda,t[i]) + Delta4(lambda,n,t[i]) ) ) - kappa_02 )
      
    }
    
  }  
  
  return(H)
  
}


D1_SingP = function(lambda, n, t, Expression =2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  
  t1 = log(q1/p1)/(q1 + p1)
  
  K_plus  = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (eta1^r)*( (2-(p1*r)/eta1) - (1-(p1*r)/eta1)*(q1/p1) )*((q1/p1)^(b_j[1]/eta1))
  K_minus = ((q1/p1)^(bj0/eta1)) - (2 - (q1/p1))*((q1/p1)^(b_j[1]/eta1))
  C0_u    = K_plus - Delta3(lambda,n,t1)
  C0_d    = K_minus- Delta1(lambda,n,t1)
  
  tilde_K_plus  = ((2*q1)^r)*((p1^2)/q1) - (eta1^r)*( p1*(2-(p1*r)/eta1) - q1*(1-(p1*r)/eta1) )
  tilde_K_minus = ((p1^2)/q1) + q1 - 2*p1
  
  H = NULL
  
  for (i in 1:L_t){
	  
    const1 = p1*((p1/q1)^(b_j[1]/eta1))*exp(-b_j[1]*t[i])
    
    if(Expression == 2){ 
      
      tilde_d1_plus  = (exp(-eta1*t[i])/bj0) + (exp(eta1*t[i])/b_j[2]) -    (2/b_j[1])
      tilde_d1_minus = (eta1^r)*( ((2 - (p1*r)/eta1 )/b_j[1]) - (( (1 - (p1*r)/eta1 )*exp(eta1*t[i]) )/b_j[2]) )
      
      kappa_0 = as.numeric(  tilde_K_plus*q2*(1-w)*( (1/bj0)  - ((b_j[2]+eta1)/(b_j[1]*b_j[2])) )*exp(-b_j[1]*t[i])  +
                               tilde_K_minus*(eta1^r)*( 1 + q2*(1-w)*((b_j[2]+eta1 - p1*r)/(b_j[1]*b_j[2])) )*exp(-b_j[1]*t[i]) ) +
        as.numeric(  const1*(C0_d*(Delta3(lambda,n,0) + q2*(1-w)*Delta4(lambda,n,0) )   - C0_u*( Delta1(lambda,n,0) + q2*(1-w)*Delta2(lambda,n,0)) -
                               ( 1 - S_1(lambda,0) + q2*(1-w)*( (1/bj0)  - S_2(lambda,0)) )*Delta3(lambda,n,t1)  - ( S_3(lambda,0) + q2*(1-w)*S_4(lambda,0)  )*Delta1(lambda,n,t1) ))   
      
      
      H = c(H, as.numeric( q2*(1-w)*(  tilde_K_plus*tilde_d1_plus +    tilde_K_minus*tilde_d1_minus ) + 
                             q2*(1-w)*const1*( C0_d*Delta4(lambda,n,t[i]) - C0_u*Delta2(lambda,n,t[i]) - ( (exp(bj0*t[i])/bj0) - S_2(lambda,t[i]))*Delta3(lambda,n,t1) - 
                                                 S_4(lambda,t[i])*Delta1(lambda,n,t1) )   ) -  kappa_0 )
    }	
    
  }  
  
  return(H)
  
}




#----------------------------------------------------------------------------------------
# Approximation of the numerator in the true function after t_max
#----------------------------------------------------------------------------------------

# Non normalized function

hat_u2_SingP = function(lambda,n,t,Expression = 2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  b0 = bj0^(0:(n-1))
  g_j0 = (q1/p1)^(b_j/eta1)
  
  inv_B = H_inv_B_SingP(lambda,n)
  u0 = der_u0_SingP(lambda,n)
  t1 = log(q1/p1)/(q1 + p1)
  
  C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (t(g_j0)%*%(inv_B%*%u0))
  C0_d = ((q1/p1)^(bj0/eta1)) - (t(g_j0)%*%(inv_B%*%b0))
  C_0  = C0_u/C0_d
  
  H = NULL
  
  for (i in 1:L_t){
    		
    if(Expression == 2){ # Delta decomposition from the expression  in the proof of Proposition 2
      
      C0_u = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - S_3(lambda,t1) - Delta3(lambda,n,t1)
      C0_d = ((q1/p1)^(bj0/eta1)) - S_1(lambda,t1) - Delta1(lambda,n,t1)
      C_0  = C0_u/C0_d
      
      tilde_b = C_0*( bj0*((q1/p1)^(bj0/eta1)) - S_5(lambda,t1) - Delta5(lambda,n,t1) ) +  S_6(lambda,t1) + Delta6(lambda,n,t1)
      
      H = c(H, as.numeric( ((2*p1)^r)*((p1/q1)^(eta2/eta1)) + tilde_b*(t[i]-t1)  ))
      
    }
    
  }  
  
  return(H)
}


N2_SingP = function(lambda,n,t, Expression = 2){ # second normalization function (multiplication by (p1*((p1/q1)^(b_j[1]/eta1)))^4  )
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  
  b_j = eta1*(1:n) + bj0
  
  t1 = log(q1/p1)/(q1 + p1)
  
  H = NULL
  
  for (i in 1:L_t){
    K_plus  = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (eta1^r)*( (2-(p1*r)/eta1) - (1-(p1*r)/eta1)*(q1/p1) )*((q1/p1)^(b_j[1]/eta1))
    K_minus = ((q1/p1)^(bj0/eta1)) - (2 - (q1/p1))*((q1/p1)^(b_j[1]/eta1))
    C0_u    = K_plus - Delta3(lambda,n,t1)
    C0_d    = K_minus- Delta1(lambda,n,t1)
    tilde_K_plus  = ((2*q1)^r)*((p1^2)/q1) - (eta1^r)*( p1*(2-(p1*r)/eta1) - q1*(1-(p1*r)/eta1) )
    tilde_K_minus = ((p1^2)/q1) + q1 - 2*p1
    const2        = (p1*((p1/q1)^(b_j[1]/eta1)))^2 
    
    if(Expression == 2){ # Delta decomposition from the expression of hat_u1 in the proof of Proposition 2
      
      H = c(H, as.numeric( tilde_K_plus*(bj0*((p1^2)/q1) - (2*p1*b_j[1]-q1*b_j[2] ))*(t[i]-t1) +  
                                     tilde_K_minus*(((2*q1)^r)*((p1^2)/q1)+ (eta1^r)*( p1*b_j[1]*(2-(p1*r)/eta1) - q1*b_j[2]*(1-(p1*r)/eta1))*(t[i]-t1) )+
                             const2*( (C0_d*Delta6(lambda,n,t1)- C0_u*Delta5(lambda,n,t1)  )*(t[i]-t1) - Delta3(lambda,n,t1) *( bj0*((q1/p1)^(bj0/eta1)) - S_5(lambda,t1) )*(t[i]-t1) -
                                     Delta1(lambda,n,t1) *( ((2*p1)^r)*((p1/q1)^(eta2/eta1)) + S_6(lambda,t1)*(t[i]-t1) ) ) ) )
      
    }
    
  }  
  
  return(H)
  
}


#----------------------------------------------------------------------------------------
# Approximation of the denominator in the true function after t_max
#----------------------------------------------------------------------------------------

D2_SingP = function(lambda,n,t,Expression = 2){
  
  L_t = length(t)
  p1 = lambda[1]
  q1 = lambda[2]
  p2 = lambda[3]
  q2 = lambda[4]
  w = lambda[5]
  w = max(1.0e-15, w)
  
  eta2 = p2+q2
  eta1 = p1+q1
  
  sicma = (eta2*q1)/(q2*w)
  r = (-q2*w)/q1
  
  bj0 = (sicma - eta1)*r
  b_j = eta1*(1:n) + bj0
  t1 = log(q1/p1)/(q1 + p1)
  
  K_plus  = ((2*p1)^r)*((p1/q1)^(eta2/eta1)) - (eta1^r)*( (2-(p1*r)/eta1) - (1-(p1*r)/eta1)*(q1/p1) )*((q1/p1)^(b_j[1]/eta1))
  K_minus = ((q1/p1)^(bj0/eta1)) - (2 - (q1/p1))*((q1/p1)^(b_j[1]/eta1))
  C0_u    = K_plus - Delta3(lambda,n,t1)
  C0_d    = K_minus- Delta1(lambda,n,t1)
  tilde_K_plus  = ((2*q1)^r)*((p1^2)/q1) - (eta1^r)*( p1*(2-(p1*r)/eta1) - q1*(1-(p1*r)/eta1) )
  tilde_K_minus = ((p1^2)/q1) + q1 - 2*p1
  const2        = (p1*((p1/q1)^(b_j[1]/eta1)))^2 
  
  H = NULL
  
  for (i in 1:L_t){

    if(Expression == 2){ 
  
      H = c(H, as.numeric( 0.5*q2*(1-w)*( tilde_K_plus*( bj0*((p1^2)/q1) - (2*p1*b_j[1]-q1*b_j[2] ) )*(t[i]-t1) + 
      tilde_K_minus*( 2*((2*q1)^r)*((p1^2)/q1) + (eta1^r)*( p1*b_j[1]*(2-(p1*r)/eta1) - q1*b_j[2]*(1-(p1*r)/eta1))*(t[i]-t1) )
      )*(t[i]-t1)  +0.5*q2*(1-w)*const2*( (C0_d*Delta6(lambda,n,t1)- C0_u*Delta5(lambda,n,t1))*(t[i]-t1) -
      Delta3(lambda,n,t1)*( bj0*((q1/p1)^(bj0/eta1)) - S_5(lambda,t1) )*(t[i]-t1) - Delta1(lambda,n,t1)*( 2*(((2*p1)^r)*((p1/q1)^(eta2/eta1))) + S_6(lambda,t1)*(t[i]-t1) )
      )*(t[i]-t1) + p1*D1_SingP(lambda,n,t1,Expression = 1)   ))
      
    }
  }  
  
  return(H)
}




#----------------------------------------------------------------------------------------
# The overall approximation function of u()
#----------------------------------------------------------------------------------------

hat_u_SingP = function(lambda,n,t, Expression = 2){  
  
	#lambda = [alpha,p1,q1,p2,q2,w]

	p1 = lambda[1]
	q1 = lambda[2]
	eta1 = p1+q1
  
	t1 = log(q1/p1)/eta1
  
	tt = t[which(t < t1)]
	f = NULL
	if(length(tt)!=0){
      f = c(f, hat_u1_SingP(lambda,n,tt, Expression= Expression ))
	}
  
	if(max(t)>= t1){
    
		ttt = t[which(t>=t1)]
		f2 = hat_u2_SingP(lambda,n,ttt, Expression= Expression )
    
		f = c(f,f2)
  }
  
  return(f)
  
}


#----------------------------------------------------------------------------------------
# The overall approximation function of F_2
#----------------------------------------------------------------------------------------

F_hat_1_SingP = function(lambda,n,t, Expression = 7){  
  
    f1 = 1 + (N1_SingP(lambda,n,t, Expression = 1)/D1_SingP(lambda,n,t, Expression = 1))
  
	return(f1)
}

#------------------------------------------------------------

F_hat_2_SingP = function(lambda, n, t, Expression = 7, delta0 = 0, delta1 = 1){  
  
	p1 = lambda[1]
	q1 = lambda[2]
	t1 = log(q1/p1)/(q1 + p1)
  
	f2 = 1 + (N2_SingP(lambda,n,t, Expression = 1)/D2_SingP(lambda,n,t, Expression = 1))

	f2 = f2 - delta0*(t-t1)*(1-delta1^(t-t1))
  
	return(f2)
}

#--------------------------------------------------------

F_hat_SingP = function(lambda,n,t, Expression =7, delta0 = 0, delta1 = 1){  
  
  p1 = lambda[1]
  q1 = lambda[2]
  eta1 = p1+q1
  t1 = log(q1/p1)/(q1 + p1)
  
  tt = t[which(t < t1)]
  
  f = c()
  if(length(tt)!=0){
      f = c(f, F_hat_1_SingP(lambda, n, tt, Expression= Expression ))
  }
  if(max(t)>=t1){
    
    ttt = t[which(t >= t1)]
    f2 = F_hat_2_SingP(lambda,n,ttt,Expression= Expression, delta0, delta1)
    
    f = c(f,f2)
  }
  
  return(f)
}

#--------------------------------------------------------

F_influential <-function(p1,q1,t){
  eta1 = p1+q1
  ( 1-exp(-rep(eta1,length(t))*t))/(1+(q1/p1)*exp(-rep(eta1,length(t))*t))
}

