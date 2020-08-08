******************************************************************************** 
** RDPOWER Stata Package 
** Empirical Illustration
** Authors: Matias D. Cattaneo, Rocio Titiunik and Gonzalo Vazquez-Bare
** Date: 03-Apr-2020
********************************************************************************
** net install rdpower, from(https://sites.google.com/site/rdpackages/rdpower/stata) replace
********************************************************************************
** NOTE: If you are using rdrobust version 2020 or newer, the option 
** "masspoints(off) stdvars(on)" may be needed in order to replicate the 
** results in the paper. For example, line 31:
**
**     rdpow demvoteshfor2 demmv, tau(5)
**
** should be replaced by:
**
**     rdpow demvoteshfor2 demmv, tau(5) masspoints(off) stdvars(on)
********************************************************************************

use rdpower_senate.dta, clear
sum demmv demvoteshfor2 population dopen dmidterm

********************************************************************************
** RDPOWER
********************************************************************************

********************************************************************************
** rdpower against tau = 5
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5)

********************************************************************************
** rdpower with covariates
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5) covs(population dopen dmidterm)

********************************************************************************
** rdpower with user-specified plot options
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5) plot graph_range(-9 9) graph_step(2) ///
		                     graph_options(title(Power function) ///
							 xline(0, lcolor(black) lpattern(dash)) ///
							 yline(.05, lpattern(shortdash) lcolor(black)) ///
							 xtitle(tau) ytitle(power) ///
							 graphregion(fcolor(white))) 

********************************************************************************
** rdpower with rdrobust options
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5) h(16 18) b(18 20)

rdpow demvoteshfor2 demmv, tau(5) kernel(uniform) vce(cluster state)

rdpow demvoteshfor2 demmv, tau(5) bwselect(certwo) vce(hc3) scaleregul(0) rho(1)

********************************************************************************
** rdpower with conventional inference
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5) all

********************************************************************************
** rdpower with user-specified bias and variance
********************************************************************************

qui rdrobust demvoteshfor2 demmv

local samph = e(h_l)
local sampsi_l = e(N_h_l)
local sampsi_r = e(N_h_r)

local bias_l = e(bias_l)/e(h_l)
local bias_r = e(bias_r)/e(h_r)

mat VL_RB = e(V_rb_l)
mat VR_RB = e(V_rb_r)

local Vl_rb = e(N)*e(h_l)*VL_RB[1,1]
local Vr_rb = e(N)*e(h_r)*VR_RB[1,1]

rdpow demvoteshfor2 demmv, tau(5) bias(`bias_l' `bias_r') ///
                             var(`Vl_rb' `Vr_rb') ///
							 samph(`samph') sampsi(`sampsi_l' `sampsi_r')

********************************************************************************
** rdpower manually increasing variance by 20%
********************************************************************************

qui rdrobust demvoteshfor2 demmv

mat VL_RB = e(V_rb_l)
mat VR_RB = e(V_rb_r)

local Vl_rb = e(N)*e(h_l)*VL_RB[1,1]*1.2
local Vr_rb = e(N)*e(h_r)*VR_RB[1,1]*1.2

rdpow demvoteshfor2 demmv, tau(5) var(`Vl_rb' `Vr_rb')

********************************************************************************
** rdpower without data
********************************************************************************

qui rdpow demvoteshfor2 demmv, tau(5)
rdpow, tau(5) nsamples(r(N_l) r(N_h_l) r(N_r) r(N_h_r)) ///
		 bias(r(bias_l) r(bias_r)) ///
		 var(r(Vl_rb) r(Vr_rb)) sampsi(r(sampsi_l) r(sampsi_r)) ///
		 samph(r(samph_l) r(samph_r))

********************************************************************************
** comparing ex-post power across specifications
********************************************************************************

rdpow demvoteshfor2 demmv, tau(5) p(1) h(20) plot

rdpow demvoteshfor2 demmv, tau(5) p(2) h(20) plot

rdpow demvoteshfor2 demmv, tau(5) p(1) plot

rdpow demvoteshfor2 demmv, tau(5) p(2) plot

********************************************************************************
** RDSAMPSI
********************************************************************************

********************************************************************************
** rsampsi with tau = 5
********************************************************************************

rdsampsi demvoteshfor2 demmv, tau(5)

********************************************************************************
** rsampsi with tau = 5 setting bandwdith and nratio with plot
********************************************************************************

rdsampsi demvoteshfor2 demmv, tau(5) beta(.9) samph(18 19) nratio(.5) plot

********************************************************************************
** rsampsi with conventional inference
********************************************************************************

rdsampsi demvoteshfor2 demmv, tau(5) all

********************************************************************************
** rsampsi vs rdpower
********************************************************************************

qui rdsampsi demvoteshfor2 demmv, tau(5)
rdpow demvoteshfor2 demmv, tau(5) sampsi(r(sampsi_h_l) r(sampsi_h_r))

********************************************************************************
** rsampsi without data
********************************************************************************

qui rdsampsi demvoteshfor2 demmv, tau(5)
local init = r(init_cond)
rdsampsi, tau(5) nsamples(r(N_l) r(N_h_l) r(N_r) r(N_h_r)) ///
				 bias(r(bias_l) r(bias_r)) ///
				 var(r(var_l) r(var_r)) ///
				 samph(r(samph_l) r(samph_r)) ///
				 init_cond(`init')

********************************************************************************
** comparing sample sizes across designs
********************************************************************************

rdsampsi demvoteshfor2 demmv, tau(5) p(0) h(20) plot

rdsampsi demvoteshfor2 demmv, tau(5) p(1) h(20) plot

rdsampsi demvoteshfor2 demmv, tau(5) p(0) plot

rdsampsi demvoteshfor2 demmv, tau(5) p(1) plot
