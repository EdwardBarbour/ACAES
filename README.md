# ACAES
Matlab scripts for running ACAES simulations
FROM PAPER  Adiabatic Compressed Air Energy Storage with packed bed thermal energy storage
E Barbour et al.
Detailed derivations can be found in the above paper.

The program simulates increments of air being compressed and added to a fixed volume store through 
packed bed heat exchangers. Each step in the simulation adds an increment of air dn to the air store, 
through the compressors and heat exchangers. The mass of air in the store increases as air as added
to the store, AS DOES THE MASS OF AIR IN EACH HEAT EXCHANGER.

During the expansion phase the reverse process is simulated, as increments of air are removed from
the high pressure air store, through the hot packed beds and then expanders (turbines). Pressure
drops are calculated using the ergun eqn. Temperature losses occur throughout the simulation from 
the packed beds.  

Instructions
Download all and run file Single_cycle_new.m

In Single_cycle_new.m 
lines 12-24 are the simulation constants
line 26 don't change the r value
lines 28-31 choose the work to be stored, number of stages, initial and final pressures
line 33 pressure drop for guessing teh final r value, guess according to previous calculated p drop
line 41-51 the packed bed heat exchanger characteristics

Comments are included in all the scripts