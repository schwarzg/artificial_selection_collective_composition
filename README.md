Artificial selection of microbial collectives with certian composition
===
## Artificial selection in simple population dynamics model 

Ref: in preparation  
This repository simulate artificial selection of microbial collectives with certain composition.  
Two-type population dynamics is assumed.  

### Simple population dynamics model
Two genotypes are present. One is considered as a mutant of the others.
Mutant grows faster.
Wildtype changed into mutant, but reversal process is not occured.

### Artificial selection protocol 
Maturation: grow population following spd model.  
Selection: Choose one collective whose mutant frequency is closest to the target value.  
Reproduction: Binomial sampling of N0 cells based on the selected collective. 

## Main codes
1. `AGS_data_generator_sampling.py`
	+ Generate trajectory of species abandance with artificial selection using sampling method.  
	+ Output data folder: `"data/raw/AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ format: array[ncycle, 1(time)+1(index)+ncomm(Initial w, Final w)+ncomm(Initial m, final m)]  
	+ Requirements : `gillespie/model.py, growth_alter.py, selection.py, mpi4py` package if possible  
1. `AGS_data_generator_stochastic.py`
	+ Generate trajectory of species abandance with artificial selection using stochastic simulation (tau-leaping algorithm).  
	+ Output data file: `"AGS_PD_sto_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ format: array[ncycle, 1(time)+1(index)+ncomm(Initial w, Final w)+ncomm(Initial m, final m)]  
	+ Requirements : `gillespie/model.py, gillespie/tau_leaping.py, selection.py, mpi4py` package if possible 
1. `NS_data_generator_stochastic.py`
	+ Generate trajectory of species abandance without artificial selection using stochastic simulation (tau-leaping algorithm).  
	+ Output data file: `"NS_PD_sto_N0%s_mbar%s_r%s_s%s_mu%s_ncycle%d"%(N0,mbar,r,s,mu,ncycle)`
	+ format: array[ncycle, 1(time)+2(Initial w, Final w)+2(Initial m, final m)]  
	+ Requirements : `gillespie/model.py, gillespie/tau_leaping.py, mpi4py` package if possible 
1. `conditional_probablity_data_generator.py`
	+ Generate conditional probability of the selected frequency using sampling simuation and (semi-)analytic expression.
	+ Output data file: `"data/one/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)`
	+ format:[n\_fk,30(bins)] 
	+ Requirements : `gillespie/model.py, analytic_results.py, growth_alter.py, selection.py, stats_custom_tool.py`

## Relative Functions
1. `gillespie`
1. `analytic_results.py`
1. `growth_alter.py`
1. `selection.py`
1. `stats_custom_tool.py`


