# elastgen 
Original author: Ralph Minehart (University of Virginia)
Modifications by K. Joo and L.C. Smith (UVA)

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

## input arguments
      call elast_gen(file1,     !LUND event file
     1               file2,     !Summary
     1               file3,     !ntuple event file
     1               file4,     !Summary
     1                   f,     !fraction of events in elas peak
     1                   d,     !egam < delta considered elastic
     1                  g_,     !Number of particles to detect / event
     1               tl,tr,     !target length, radius (cm)
     1             vx,vy,vz,    !target vertex coordinates
     1                  be,     !beam energy
     1                 emn,     !minimum scattered e- energy
     1             smn,smx,     !min,max e- scattering angle
     1                  p_,     !0=include elastic 1=exclude
     1                   M,     !Process M events
     1                   c,     !Maximum cross section
     1  opt_strag,opt_fiduc,    !Straggling of scattered e- (1=on,0=off) and fiducial cut options
     1             idum_off)    !random number offset

## demo input file (inp/elas_clas12.inp)
     0.2
     0.005
     3
     4.0 0.6
     0.0 0.0 0.0 
     10.600
     2.5
     5. 20.
     0
     100000
     0.00002
     1,0
     123

