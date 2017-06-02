
##Examples
library(QUIC)
set.seed(1)
x1 = matrix(rnorm(50 * 20),ncol = 20)
x2 = matrix(rnorm(50 * 20),ncol = 20)
z = matrix(rnorm(50 * 20),ncol = 20)
Y = array(0, c(2, 50, 20))
Y[1, , ] = x1 + z
Y[2, , ] = x2 + z

a = EMgCov(Y, lambda.i = 0.1, lambda.z = 0.02)

EMgCov <- function(Y, lambda.i = 0.4, lambda.z = 0.1, 
		 Omega.threshold = 0.005, max_iter = 100, tol_value = 0.01, eigen.threshold = 0.01,
		 alpha.value = rep(1, 4), Use.alpha = FALSE,  Use.QUIC = TRUE, Inf.Projection = FALSE,
		 pen.diag = F) {
	#########################################
	#### Arguments ##########################
	#########################################
	# This is estimate the dependent graph using EM algorithm
	# If Use.alpha = TRUE, then use the extended version Y_k = X_k + alpha Z
	# Y is a K x n x p array where K is the number of categories, n is sample size, p is dimension 
	# lambda.i is the tuning parameter for categorical networks
	# lambda.z is the tuning paraemter for the systemic network
	# Omega.threshold is the hard thresholdind parameter for the final results.
	# max_iter is the number of maximum iterations
	# tol_value control the convergence rate
	# eigen.threhold is the parameter controlling the eigen values for the intial Omega's
	# alpha.value is the tuning parameter for the extended version
	# Inf.Projection is control the method used to intialize the algorithm
	#########################################
	######## Value ##########################
	#########################################
	# Omega is a K x p x p array corresponding to the estimated precision matrices using EM 
	# Omega.One is a K x p x p array corresponding to the estimated precision matrices using One step method
	# alpha.value is the estimated alpha value for the extended model
	count = 1
	diff_value = 1
	K = dim(Y)[1]
	n = dim(Y)[2]
	p = dim(Y)[3]
	e = new.env(hash = T)
	rho.matrix = matrix(1, p, p)
	diag(rho.matrix) = 0
	
	## Set the optimizaiton parameters  
	
	Omega_new = array(0, c((K + 1), p, p))
	e$Sigma = array(0, c((K + 1), p, p))
	
	Omega = GetOmega(Y = Y, lambda.i = lambda.i, lambda.z = lambda.z, 
			eigen.threshold = eigen.threshold, Inf.Projection = Inf.Projection,
			alpha.value = alpha.value, Use.QUIC = Use.QUIC)

	Omega.One = Omega
	cat(sep = "", "[", 0,"]","dif=", round(diff_value, 3), "\n")
	##################
	##start loop
	##################
	while((diff_value > tol_value) & (count < max_iter))
	{	
		A = Omega[K + 1, , ]
		for(k in 1:K){
			A = A + alpha.value[k]^2 * Omega[k, , ]
		}
		A0 = A - Omega[K + 1, , ]
		A.Chol = solve(chol(A))
		A.inv = A.Chol %*% t(A.Chol)
		# B is sum of Omega_k Sigma_k Omega_k
		C = matrix(0, p, n)
		
		for(k in seq(1, K)){
			C = C + alpha.value[k] * Omega[k, , ] %*% t(Y[k, , ])
		}
		
		Compute.Sigma(e, A.inv = A.inv, C = C, 
				Y = Y, Omega = Omega, A0 = A0, alpha.value = alpha.value)
		
		for(k in seq(1, K)){ 
			if(Use.QUIC){
				Omega_new[k, , ] = QUIC(e$Sigma[k, , ], rho = lambda.i * rho.matrix)$X 
			} else {
				Omega_new[k, , ] = glasso(e$Sigma[k, , ], rho = lambda.i, thr = Omega.threshold,
						maxit = 100, penalize.diagonal= pen.diag)$wi   
				Omega_new[k, , ] = (Omega_new[k, , ] + t(Omega_new[k, , ]))/2	        
			}
		}
		
		if(Use.QUIC){
			Omega_new[K + 1, , ] = QUIC(e$Sigma[K + 1, , ], rho = lambda.z * rho.matrix)$X       
		}else {
			Omega_new[K + 1, , ] = glasso(e$Sigma[K + 1, , ], rho = lambda.z, thr = Omega.threshold,
					maxit = 100, penalize.diagonal= pen.diag)$wi  
			Omega_new[K + 1, , ] = (Omega_new[K + 1, , ] + t(Omega_new[K + 1, , ]))/2      
		}
		
		if(Use.alpha) {
			for(k in 1:K){
				Ez = A.inv %*% C 
				Ezz = n * A.inv + Ez %*% t(Ez) 
				Y_k = Y[k, , ]
				Sk = (t(Y_k) %*% t(Ez) +  Ez %*% Y_k) /(2 * n)
				alpha.value[k] = tr(as.matrix(Omega_new[k, , ] %*% Sk)) / tr(as.matrix(Omega_new[k, , ] 
										%*% Ezz /n ))
			}
		}
		
		# diif_value should be compute before assigning OMEGA
		diff_value = sum(abs(Omega[ , , ] - Omega_new[ , , ]))/ sum(apply(abs(Omega), 1, sum))
		Omega = Omega_new
		
		cat(sep = "", "[",count,"]","dif=", round(diff_value, 3), "\n")
		
		if(diff_value < tol_value)
		{
			cat(sep = "" ,"\n")
		}
		count = count + 1
	}
	Thresh.id = (abs(Omega) < Omega.threshold)
	Omega[Thresh.id] = 0
	if(Use.alpha){
		Omega.scale =max(diag(solve(Omega[K + 1, , ])))
		Omega[K + 1, , ] =  Omega[K + 1, , ] * Omega.scale
		alpha.value = alpha.value * Omega.scale
	}
	
	result = list(Omega = Omega, Omega.One = Omega.One, alpha.value = alpha.value)
	gc()
	return(result)
}

GetOmega <- function(Y, lambda.i=0.4, lambda.z=0.1, pen.diag = FALSE,
		eigen.threshold = 0.01, alpha.value = rep(1, 4), Use.QUIC = FALSE,
		Inf.Projection = FALSE){
	K = dim(Y)[1]
	n = dim(Y)[2]
	p = dim(Y)[3]
	e = new.env(hash = T)
	rho.matrix = matrix(1, p, p)
	diag(rho.matrix) = 0
	Omega <- array(0, c(K + 1, p, p))
	rho.matrix = matrix( rep(1, p^2), p, p)
	diag(rho.matrix) = 0
	
	S = GetSigma(Y, eigen.threshold, alpha.value = alpha.value, 
			Inf.Projection = Inf.Projection)       
	for(k in 1:K){
		if(Use.QUIC){
			Omega[k, , ] = QUIC(S[k, , ], rho = lambda.i * rho.matrix)$X   # 33.511 second
		}else {
			a = glasso(S[k, , ], rho = lambda.i, penalize.diagonal = pen.diag, maxit = 100)
			Omega[k, , ] = (a$wi + t(a$wi) )/2         
		}
	}
	
	if(Use.QUIC){
		Omega[K + 1, , ] = QUIC(S[K + 1, , ], rho = lambda.z * rho.matrix)$X  
	}else {
		a = glasso(S[K + 1, , ], rho = lambda.z, penalize.diagonal = pen.diag, maxit = 100)
		Omega[(K + 1), , ] = (a$wi + t(a$wi) )/2
	}
	return(Omega)
}


Compute.Sigma <- function(e, A.inv, C, Y, Omega, A0, alpha.value = rep(1,4)) 
{
	K = dim(Y)[1]
	n = dim(Y)[2]
	p = dim(Y)[3]
	B = matrix(0, p, p) 
	Ez = A.inv %*% C 
	Ezz = n * A.inv + Ez %*% t(Ez)   ##### be careful here for the n
	for(k in 1:K){
		EYZk = alpha.value[k] * t(Y[k, , ]) %*% t(Ez)
		e$Sigma[k, , ] = as.matrix((t(Y[k, , ]) %*% Y[k, , ] - EYZk - t(EYZk) +
							alpha.value[k]^2 * Ezz) / n) 
	}
	e$Sigma[K + 1, , ] = as.matrix(Ezz / n)
}

GetSigma = function(Y, eigen.threshold = 0.01, alpha.value = rep(1, 4), 
		Inf.Projection = FALSE){
	K = dim(Y)[1]
	n = dim(Y)[2]
	p = dim(Y)[3]
	S = array(0, c(K + 1, p, p))
	temp = matrix(0, p, p)   
	
	for(i in 1:(K-1)){
		for(j in (i+1):K){
			temp = temp + cov(Y[i, , ], Y[j, , ])
		}
	}
	
	alpha.value = as.matrix(alpha.value)
	alpha.weigh = sum( alpha.value %*% t(alpha.value)) -  t(alpha.value) %*% alpha.value 
	temp = (temp + t(temp)) 
	temp = temp / alpha.weigh[1]
	S[K + 1, , ] = temp
	
	if(eigen(S[K + 1, , ])$values[p] < 0){
		if(Inf.Projection){
			S[K + 1, , ] = MatrixMaxProj(S[K + 1, , ], tol_value = 0.01)
		}else{
			S[K + 1, , ] = GetSigX(S[K + 1, , ], eigen.threshold)
		}
	}
	
	for(k in 1:K){
		S[k, , ] = var(Y[k, , ]) - alpha.value[k, 1]^2 * S[K + 1, , ]
		if(eigen(S[k, , ])$values[p] < 0){
			if(Inf.Projection){
				S[k, , ] = MatrixMaxProj(S[k, , ], tol_value = 0.01)
			}else{
				S[k, , ] = GetSigX(S[k, , ], eigen.threshold)
			}
		}
	}    
	return(S)
}

GetSigX = function(A, eigen.threshold = 0.01){
	# this function is to make sure the matrix have lowest positive eigen value
	# This will keep the same eigen vector but only threshold the eigen value to
	# a presepecified eigen value
	eig.A = eigen(A)$value
	id = (eig.A < eigen.threshold)
	eig.A[id] = eigen.threshold
	result = eigen(A)$vectors %*% diag(eig.A) %*% t(eigen(A)$vectors)
	return(result)	
}

MatrixMaxProj = function(M, r = 0.5, tol_value = 0.01, N = 100) {
	#This is the projection to make the matrix positive definite
	t = 0  
	p = dim(M)[1]
	Rt = GetSigX(M, 0.01)
	Zt = matrix(rep(1/p^2, p^2), p, p)
	while(t < N){ 
		R0t = Rt - GetSigX(Rt - Zt, eigen.threshold = 0) # ex
		Z0t = Zt - PB1.projection(Rt + Zt - M)     # ez
		diff = max(max(abs(R0t)), max(abs(Z0t)))
		cat(sep = "", "[",t,"]","dif=", round(diff, 3), "\n")
		if(diff < tol_value){
			break
		}
		Rt = Rt - r * (R0t - Z0t)/2
		Zt = Zt - r * (R0t + Z0t)/2
		t = t + 1
	}
	A = as.matrix(Rt)
	return(A)
}

PB1.projection = function(A){
	# This is to find matrix A.result where |A.result|_1 <=1
	# min ||  A - A.result||_F
	p = dim(A)[1]
	a = as.vector(A)
	T.mat = T.matrix(A * sign(A))
	a1 = T.mat %*% (sign(a) * a)  
	a.del = Delta.vec(a1)
	a.del.cum = Cumulate.vec( a.del * (1:p^2) )
	b = Cumulate.vec(a1)
	if(b[p^2] <= 1){
		x = a1    
	} else {
		id = (a.del.cum < 1)
		K = sum(id) + 1
		y = (b[K] - 1)/K
		x = rep(0, p^2)
		x[1:K] = a1[1:K] - y
	}
	B1 = sign(a) * solve(T.mat) %*% x
	A.result = matrix(B1, p, p)
	A.result = (A.result + t(A.result))/2
	return(A.result)
}

T.matrix = function(A){
	p = dim(A)[1]
	a = as.vector(A)  
	id = order(a, decreasing = TRUE)
	T.mat = matrix(0, p^2, p^2)
	mat.id = cbind(1:p^2, id)
	T.mat[mat.id] = 1  
	return(T.mat)
}

Delta.vec = function(a){
	# calculate the different between consective entries
	p = length(a)
	a.result = rep(0, p)
	a.result[p] = a[p]
	a.result[1:(p - 1)] = a[1:(p - 1)] - a[2:p]
	return(a.result)
}

Cumulate.vec = function(a){
	# calculate the cumulate sum for a vector
	p = length(a)
	b = rep(0, p)
	for(i in 1:p){
		b[i] = sum(a[1:i])
	}
	return(b)
}

