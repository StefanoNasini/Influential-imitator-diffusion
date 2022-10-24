
rm(list=ls(all=TRUE))
gc(reset = TRUE)

library(Rmpfr)

###############################################################################################
###############################################################################################
# 
# Load functions
#
###############################################################################################
###############################################################################################

source("DiffusionFun.R")

#----------------------------------------------------------
# Generate plots
#----------------------------------------------------------

CASES = c(0.01, 0.05, 0.1, 0.2)

w = 0.5
N = 5
t = seq(0,6,by=0.1)

for(cc in CASES){ 

	p1 = cc     
	q1 = 2                         
	p2 = cc
	q2 = 2                           
	lambda=c(p1,q1,p2,q2,w)
	(t_max = log(q1/p1)/(q1 + p1))
	
	filenamePLOT = paste("Plot_U_approx", as.character(100*cc),  ".pdf", sep = "")

	pdf(filenamePLOT)
	plot(t,Num_True_SingP(lambda,t),type = 'l', xlab = "time", ylab = "u(t)")
	
	for(n in 2:N){
		lines(t,hat_u1_SingP(lambda,n,t),col= "red", lty = n, lwd = n-1)
	}

	abline(v = as.numeric(t_max), col = "blue", lwd = 3)
	dev.off()
}


#----------------------------------------------------------------------------------------

for(cc in CASES){ 

	p1 = cc   
	q1 = 1-p1                         
	p2 = cc
	q2 = 1-p2                            
	lambda=c(p1,q1,p2,q2,w)
	(t_max = log(q1/p1)/(q1 + p1))
	
	filenamePLOT = paste("Plot_F_approx", as.character(100*cc),  ".pdf", sep = "")

	pdf(filenamePLOT)
	plot(t,F_True_SingP(lambda,t),type = 'l', xlab = "time", ylab = "F(t)")
	
	for(n in 2:N){
		lines(t,F_hat_1_SingP(lambda,n,t),col= "red", lty = n, lwd = n-1)
	}

	abline(v = as.numeric(t_max), col = "blue", lwd = 3)
	dev.off()
}


#----------------------------------------------------------------------------------------


for(cc in CASES){ 

	p1 = cc   
	q1 = 1-p1                         
	p2 = cc
	q2 = 1-p2                            
	lambda=c(p1,q1,p2,q2,w)
	(t_max = log(q1/p1)/(q1 + p1))
	
	filenamePLOT = paste("Plot_U_approx_complete", as.character(100*cc),  ".pdf", sep = "")

	pdf(filenamePLOT)
	plot(t,Num_True_SingP(lambda,t),type = 'l', xlab = "time", ylab = "u(t)")
	
	for(n in 2:N){
		lines(t,hat_u_SingP(lambda,n,t),col= "red", lty = n, lwd = n-1)
	}

	abline(v = as.numeric(t_max), col = "blue", lwd = 3)
	dev.off()
}




#----------------------------------------------------------------------------------------


for(cc in CASES){ 

	p1 = cc   
	q1 = 1 - p1                         
	p2 = cc
	q2 = 1 - p2                            
	lambda=c(p1,q1,p2,q2,w)
	(t_max = log(q1/p1)/(q1 + p1))
	
	filenamePLOT = paste("Plot_F_approx_complete", as.character(100*cc),  ".pdf", sep = "")

	pdf(filenamePLOT)
	plot(t,F_True_SingP(lambda,t),type = 'l', xlab = "time", ylab = "F(t)")
	
	for(n in 2:N){
		lines(t,F_hat_SingP(lambda,n,t, delta0 = 0.33, delta1 = 0.85), col= "red", lty = n, lwd = n-1)
	}

	abline(v = as.numeric(t_max), col = "blue", lwd = 3)
	dev.off()

}





