# elastgen     
This program makes event files (plain text and LUND format) for (e,e'p) elastic scattering
with the elastic peak modified by vertex corrections according to the Mo-Tsai 
formula (II.6) and internal bremmstrahlung (radiative tails) calculated from 
a Monte-Carlo sampling of formula (B.4) and (B.5) from RMP, Vol. 41, 205 (1969). 

Some tricks are employed to make the monte-carlo concentrate
on the regions where the cross section is large.  Also the
code attempts to find the maximum value of the radiative integrand
to normalize the monte-carlo sampling (sigr_max).  Note that 
sigr_max may be underestimated when incident e- straggling is turned
on.  The code attempts to compensate for this by oversampling those
events with large incident e- energy loss (mcall>1), which may cause 
spurious statistical spikes when f=1 and delta=0. For that choice of
options, input e- straggling should be commented out.  
