# ACAES
Matlab scripts for running ACAES simulations <br>
From the paper Adiabatic Compressed Air Energy Storage with packed bed thermal energy storage
E Barbour et al.
Detailed derivations can be found in the paper.

<p> The program simulates increments of air being compressed and added to a fixed volume store through 
packed bed heat exchangers. Each step in the simulation adds an increment of air dn to the air store, 
through the compressors and heat exchangers. The mass of air in the store increases as air as added
to the store, as does the mass of air in each of the packed beds. </p>

<p> During the expansion phase the reverse process is simulated, as increments of air are removed from
the high pressure air store, through the hot packed beds and then expanders (turbines). Pressure
drops are calculated using the ergun eqn. Temperature losses occur throughout the simulation from 
the packed beds.  </p>

<h2>Instructions</h2>
<p>Download all and run file Single_cycle_new.m<br>
  Comments are included in all the scripts</p>

<p>In Single_cycle_new.m:
  <ul>
    <li>lines 12-24 are the simulation constants</li>
    <li>line 26 don't change the r value</li>
    <li>lines 28-31 choose the work to be stored, number of stages, initial and final pressures</li>
<li>line 33 pressure drop for guessing the final r value, guess according to previous calculated p drop</li> 
    <li>line 41-51 the packed bed heat exchanger characteristics</li> 
    </ul>
  </p>

