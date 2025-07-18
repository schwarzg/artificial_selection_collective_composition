Artificial selection of microbial collectives with certian composition
===
## Artificial selection in simple population dynamics model 

Ref: https://doi.org/10.7554/eLife.97461.2

Juhee Lee, Wenying Shou, Hye Jin Park (2024) Rafting a waterfall: Artificial selection for collective composition can succeed or fail depending on the initial and target values *eLife* **13**:RP97461.

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

## Main codes - Data Generator
1. `AGS_data_generator_sampling.py`
	+ Generate trajectory of species abandance with artificial selection using sampling method.  
	+ Output data folder: `"data/raw/AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ format: array[ncycle, 1(time)+1(index)+2ncomm(Initial w, Final w)+2ncomm(Initial m, final m)]  
	+ Requirements : `gillespie/model.py, growth_alter.py, selection.py, mpi4py` package if possible  
1. `AGS_data_generator_stochastic.py`
	+ Generate trajectory of species abandance with artificial selection using stochastic simulation (tau-leaping algorithm).  
	+ Output data file: `"AGS_PD_sto_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ format: array[ncycle, 1(time)+1(index)+2ncomm(Initial w, Final w)+2ncomm(Initial m, final m)]  
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
1. `AGS_data_generator_three_stochastic.py`
	+ Generate trajectory of species abandance with artificial selection using stochastic simulation (tau-leaping algorithm).  
	+ Three genotype version
	+ Output data file: `"data/raw/AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%d_mhat%s_vhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,mhat,vhat,ncycle)`
	+ format: array[ncycle, 1(time)+1(index)+2ncomm(Initial w, Final w)+2ncomm(Initial m, final m)+2ncomm(Initial v, Final v)]  

## Main code - Data processing codes
1. `Ensemble_average.py`
	+ Generate ensemble-averaged data and error data at the end of simulation
	+ Work only for two-genotype population
	+ input data file: `"data/raw/AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ output data file: `'N0%s_r%s_s%s_mu%s_g%s_diagram_fine_dist.txt'%(N0,r,s,mu,ncomm)`,`"%s_%d_%s_%s_%s_%s_%s_AGS_nens%d.cycle_v2"%(N0,mbar,r,s,mu,ncomm,rhat,nens)`

## Plotting code
1. `plot_fig2.py`
	+ Input data folder: `"AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)` 
	+ Requirements : `gillespie/model.py, gillespie/tau_leaping.py, selection.py`
1. `plot_fig3.py` 
	+ Draw conditional probabilities and boxplot in Fig4.
	+ Input data folder: `"data/one/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)`
1. `plot_fig4.py`
	+ Draw successful region diagram in Fig4.
	+ Input data folder: `"data/one/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)`
1. `plot_fig5.py`
	+ Draw schematic figures and ternary diagram in Fig5
	+ Input data : `"data/raw/AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%s_mhat%s_vhat%s_ncycle%d%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,ncycle,e)`

1. `plot_app1fig1.py`
	+ Draw simulation results in App1 Fig1.
1. `plot_app1fig4.py`
	+ Draw heatmap of absolute error in App1 Fig4. App4 Fig1, App6 Fig1 are also drawn with this code.
	+ Input data file: `"data/ens/N0%s_r%s_s%s_mu%s_g%s_ncycle%d_diagram_abs.txt"%(N0,r,s,mu,ncomm,ncycle)`
1. `plot_app1fig3.py`
	+ Draw trajectories in App1 Fig3.
1. `plot_app2fig1.py`
	+ Draw conditional probability distribution violin plot 
	+ Input data file: `"data/one/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)`
1. `plot_app2fig2.py`
	+ Draw proportionality of variance to Newborn size.
1. `plot_app2fig3.py`
	+ Draw distribution of Adult's frequencies and difference between mean and median.
1. `plot_app7fig1.py`
	+ Draw trajectories of top and top5% tactics in FigS8.
	+ Inpot data: `"data/ens/AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle)`
	+ Inpot data: `"data/ens/AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_nsel%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,nsel,rhat,ncycle)`
1. `plot_app6fig1.py`
	+ Draw conditional probability, and heatmap when deletereious mutation
1. `plot_app8fig1.py`
	+ Draw the schematics to get conditional extreme value distribution in three-type population, computed accessible region, and the simulation results in ternary plot
	+ Input data file: `"data/raw/AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%d_mhat%s_vhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,mhat,vhat,ncycle)`


## Relative Functions
1. `gillespie`
1. `analytic_results.py`
1. `growth_alter.py`
1. `selection.py`
1. `stats_custom_tool.py`


