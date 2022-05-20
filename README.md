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
1. `AGS_data_generator.py`
1. `NS_data_generator.py`
1. `transition_probablity_sampsim.py`

## Relative Functions
1. `gillespie`
1. `analytic_results.py`
1. `growth_alter.py`
1. `reproduction.py`
1. `selection.py`


