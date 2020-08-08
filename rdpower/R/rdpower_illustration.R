###################################################################
# rdpower: power calculations for RD designs
# Illustration file
# !version 0.5 02-Apr-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################
## NOTE: if you are using rdrobust version 2020 or newer, the option 
## masspoints="off" and stdvars="on" may be needed in order to replicate the 
## results in the paper.
## For example, line 39:
##    aux = rdpower(data=Z,tau=5)
## should be replaced by:
##    aux = rdpower(data=Z,tau=5,masspoints="off",stdvars="on")
###################################################################

#################################################################
# Setup
#################################################################

rm(list = ls())

library(rdpower)
library(rdrobust)

data = read.csv('rdpower_senate.csv')
Y = data$demvoteshfor2
R = data$demmv
names(data)
# Outcome and running variables 
Z = data[c('demvoteshfor2','demmv')]
# Covariates
covs = data[c('population','dopen','dmidterm')]


#################################################################
# rdpower against tau = 5
#################################################################

aux = rdpower(data=Z,tau=5)

#################################################################
# rdpower with covariates
#################################################################

aux = rdpower(data=Z,tau=5,covs=covs)

#################################################################
# rdpower with plot
#################################################################

aux = rdpower(data=Z,tau=5,plot=TRUE)

#################################################################
# rdpower with rdrobust options
#################################################################

aux = rdpower(data=Z,tau=5,h=c(16,18),b=c(18,20))
aux = rdpower(data=Z,tau=5,kernel='uniform',cluster=data$state,vce='hc0')
aux = rdpower(data=Z,tau=5,bwselect='certwo',vce='hc3',scaleregul=0,rho=1)

#################################################################
# rdpower with conventional inference
#################################################################

aux = rdpower(data=Z,tau=5,all=TRUE)

#################################################################
# rdpower with user-specified bias and variance
#################################################################

rdr = rdrobust(Y,R)

samph = rdr$bws[1]

sampsi.l = rdr$N_h_l
sampsi.r = rdr$N_h_r

bias.l = rdr$bias[1]
bias.r = rdr$bias[2]

VL = rdr$V_rb_l
VR = rdr$V_rb_r

N = sum(rdr$N)

Vl.rb = N*rdr$bws[1]*VL[1,1]
Vr.rb = N*rdr$bws[1]*VR[1,1]

aux = rdpower(data=Z,tau=5,bias=c(bias.l,bias.r),variance=c(Vl.rb,Vr.rb),samph=samph,sampsi=c(sampsi.l,sampsi.r))

#################################################################
# rdpower manually increasing variance by 20%
#################################################################

rdr = rdrobust(Y,R)

VL = rdr$V_rb_l
VR = rdr$V_rb_r

N = sum(rdr$N)

Vl.rb = N*rdr$bws[1]*VL[1,1]*1.2
Vr.rb = N*rdr$bws[1]*VR[1,1]*1.2

aux = rdpower(data=Z,tau=5,variance=c(Vl.rb,Vr.rb))

#################################################################
# rdpower without data
#################################################################

aux1 = rdpower(data=Z,tau=5)

aux = rdpower(tau=5,
              nsamples=c(aux1$N.l,aux1$Nh.l,aux1$Nh.r,aux1$N.r),
              bias=c(aux1$bias.l,aux1$bias.r),
              variance=c(aux1$Vl.rb,aux1$Vr.rb),
              sampsi=c(aux1$sampsi.l,aux1$sampsi.r),
              samph=c(aux1$samph.l,aux1$samph.r))

#################################################################
# comparing exp-post power across specifications
#################################################################

aux1 = rdpower(data=Z,tau=5,p=1,h=20,plot=TRUE)
aux2 = rdpower(data=Z,tau=5,p=2,h=20,plot=TRUE)
aux3 = rdpower(data=Z,tau=5,p=1,plot=TRUE)
aux4 = rdpower(data=Z,tau=5,p=2,plot=TRUE)

#################################################################
## RDSAMPSI
#################################################################

#################################################################
# rdsampsi with tau = 5
#################################################################

aux = rdsampsi(data=Z,tau=5)

#################################################################
# rdsampsi setting bandwidth and nratio with plot
#################################################################

aux = rdsampsi(data=Z,tau=5,beta=.9,samph=c(18,19),nratio=.5,plot=TRUE)

#################################################################
# rdsampsi with conventional inference
#################################################################

aux = rdsampsi(data=Z,tau=5,all=TRUE)

#################################################################
# rdsampsi vs rdpower
#################################################################

aux1 = rdsampsi(data=Z,tau=5)
aux2 = rdpower(data=Z,tau=5,sampsi=c(aux1$sampsi.h.l,aux1$sampsi.h.r))

#################################################################
# rdsampsi without data
#################################################################

aux1 = rdsampsi(data=Z,tau=5)
aux2 = rdsampsi(tau = 5,
                nsamples = c(aux1$N.l,aux1$Nh.l,aux1$Nh.r,aux1$N.r),
                bias = c(aux1$bias.l,aux1$bias.r),
                variance = c(aux1$var.l,aux1$var.r),
                samph = c(aux1$samph.l,aux1$samph.r),
                init.cond = aux1$init.cond)

#################################################################
# comparing sample sizes across designs
#################################################################

aux1 = rdsampsi(data=Z,tau=5,p=0,h=20,plot=TRUE)
aux2 = rdsampsi(data=Z,tau=5,p=1,h=20,plot=TRUE)
aux3 = rdsampsi(data=Z,tau=5,p=0,plot=TRUE)
aux4 = rdsampsi(data=Z,tau=5,p=1,plot=TRUE)
