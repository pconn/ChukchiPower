# ChukchiPower
Code to conduct a power analysis of instrumental based surveys in the Chukchi Sea

To conduct an example analysis, use the "sim_driver_Chukchi_power.R" script in the ./ChukchiPower/inst directory.
Note that thihs example analysis does not conduct the whole analysis reported in the paper as this would take weeks!  Rather,
it performs a single simulation for each species and for a single effort allocation strategy and estimation model.  
This file can be edited to (roughly) reproduce our work by varying the flights, models, and number of simulations conducted.  
I say "roughly" because the particular random number generator seeds used differed from what is included in this package (which
was necessary for batch processing).

