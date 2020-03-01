# fitting r to the data

#the solution to dN/dt = b*N*(1-N/K):

logisticModel = function(t, t_0, N_0, b, K){
 K*N_0*exp(b*(t-t_0))/ (K + N_0* (exp(b*(t-t_0))- 1) )
}

#t_0 is the time the experiment began (take t_0 = 0 below)
#N_0 is the initial number of cells at time t_0
#b is the initial growth rate of the cells
#K is the carrying capacity
#we'll fit this for each experiment (dose & celltype)

#define the objective function to be the sum of squared errors
#under this 'logisticModel'


getSSE = function(drug,dose,celltype,b,K){
 if(drug == 'mtor') dat = mtor.agg
 if(drug == 'ribo') dat = ribo.agg
 inds = which(dat$DoseNum==dose & dat$Cells==celltype)
 y_k = exp(dat[inds,]$logCellNum) #fit to number of cells instead of ln(cells)
 K = K
 b=b
 logMod = function(t) logisticModel(t,0,y_k[1],b,K)
 x_k = c(y_k[1], Vectorize(logMod)(dat[inds[-1],]$Day) )
 return( sum( (x_k-y_k)^2 ) )
}


#find b that minimizes getSSE for a given drug
# and experiment (dose and celltype)

get_r = function(drug,dose,celltype,K){
 obj = function(b) getSSE(drug,dose,celltype,b,K)  #fix drug, dose, celltype, K
 opt = optim(2,obj, control=list(maxit=2000), method='Brent',lower=0,upper=2)
 b = c(opt$par)
 if(opt$convergence != 0){print('error!')} #throw an error if optim fails
 return(b)
}

# 'logisticModel' requires an initial number of cells
# and this function pulls that from the relevant experiment

getN_0 = function(drug,dose,celltype){
 if(drug == 'mtor') dat = mtor.agg
 if(drug == 'ribo') dat = ribo.agg
 ind = which(dat$DoseNum==dose & dat$Cells==celltype)[1]
 return(exp(dat[ind,]$logCellNum))
}

# returns initial number of cells, r, and K for a given
# experiment

get_N0_and_b = function(drug,dose,celltype,K){
 N_0 = getN_0(drug,dose,celltype)
 b = get_r(drug,dose,celltype,K)
 return(c(N_0,b))
}



