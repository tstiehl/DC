# Models to study DC kinetics 
This repository contains the MATLAB scripts used to simulate DC kinetics and to estimate parameters.
 - Model 1: Random entry, exit and proliferation of DCs
 - Model 2: Exit of DCs after fixed residence time, same residence time for all DCs, proliferation rate increases during the stay in the SI
 - Model 3: Exit of DCs after a maximal residence time or before, proliferation rate increases during the stay in the SI
## Executable scripts
 - main_Model1.m: generates fit and simulation of Model 1 shown in Fig. 4C
 - main_Model2.m: generates fit and simulation of Model 2 shown in Fig. 4C
 - main_Model3.m: generates fit and simulation of Model 3 shown in Fig. 4C
 - main_Model3_loop.m: generates fits and simulations of Model 3 shown in Supplementary Fig. 4D
## Scripts containing functions
- obj_Model1.m: weighted least square cost function for fitting of Model 1
- obj_Model2.m: weighted least square cost function for fitting of Model 2
- obj_Model3.m: weighted least square cost function for fitting of Model 3
- photoconverted_Model1.m: function to calculate the fraction of photoconverted cells according to Model 1 (implementation of equation 13)
- photoconverted_Model2.m: function to calculate the fraction of photoconverted cells according to Model 2 (implementation of equation 58)
- photoconverted_Model3.m: function to calculate the fraction of photoconverted cells according to Model 3 (implementation of equation 77)
## Scripts containing data
- dataPC.m: data used for fitting
