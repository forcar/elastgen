# [elastgen](https://github.com/forcar/elastgen/blob/master/src/elast_gen.F)
Original author: Ralph Minehart (University of Virginia)
Modifications by K. Joo and L.C. Smith (UVA)

This program makes event files (plain text and LUND format) for (e,e'p) elastic scattering
with the elastic peak modified by vertex corrections according to the Mo-Tsai 
formula (II.6) and internal bremmstrahlung (radiative tails) calculated from 
a Monte-Carlo sampling of formula (B.4) and (B.5) from [RMP, Vol. 41, 205 (1969)](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf)

Some tricks are employed to make the monte-carlo sampling concentrate
on the regions where the cross section has large radiative peaks.  Also the
code attempts to find the maximum value of the radiative integrand
to normalize the monte-carlo sampling (sigr_max).  Note that 
sigr_max may be underestimated when incident e- straggling is turned
on (opt_strag_0=1).  The code attempts to compensate for this by oversampling those
events with large incident e- energy loss (mcall>1), which may cause 
spurious statistical spikes when f=1 and d=0. For that choice of
options, user should set opt_strag_0=0.  

External straggling options estimate e- energy loss only and do not generate photons.

Since GEMC does not simulate beam e- straggling preferred options is to set opt_strag_0=1
and let GEMC/GEANT handle scattered e- straggling (opt_strag_1=0).

## input arguments
      call elast_gen(file1,     !LUND event file
     1               file2,     !Summary 1 
     1               file3,     !Summary 2
     1                   f,     !fraction of events in elas peak
     1                   d,     !egam < d considered elastic
     1                  g_,     !output particles 2=e-, p+ 3=e-, p+ and gamma
     1               tl,tr,     !target length, radius (cm)
     1             vx,vy,vz,    !target vertex coordinates
     1                  be,     !beam energy (Gev)
     1                 emn,     !minimum scattered e- energy
     1             smn,smx,     !min,max e- scattering angle
     1                  p_,     !0=include elastic 1=exclude
     1                   M,     !Process M events
     1                   c,     !Maximum cross section (sigr_max in code)
     1         opt_strag_0,     !Straggling of beam e- (1=on, 0=off)
     1         opt_strag_1,     !Straggling of scattered e- (1=on, 0=off)
     1           opt_fiduc,     !use fiducial cut (1=on, 0=off)
     1            idum_off)     !random number offset

## demo input file (inp/elas_clas12.inp)
     0.2
     0.005
     3
     5.0 0.6
     0.0 0.0 0.0 
     10.600
     2.5
     5. 20.
     0
     100000
     0.00002
     1,0
     0
     123
## build, run and view (ifarm)
    ifarm1801.jlab.org> scons
    ifarm1801.jlab.org> ./elastgen < inp/elas_clas12.inp
    ifarm1801.jlab.org> cd paw
    ifarm1801.jlab.org> paw
    CLASPAW> elastgen
    CLASPAW> elastgen#plot
    
<img width="743" alt="Screen Shot 2023-01-27 at 5 28 53 PM" src="https://user-images.githubusercontent.com/10797791/215218424-e7cbeb8a-6fac-4700-9250-a78cc3da4c2b.png">

## ntuple 10  
     es: incident e- energy 
     ep: scattered e- energy
     the: scattered e- theta
     w:  true hadronic invariant mass
     ppx,ppy,ppz: components of p+ momentum
     eprot: proton energy  
     thpr: proton theta
     phpr: proton phi
     mm: missing mass
     cstk,phk: photon polar angles relative to the q vector
     eg: photon energy
     egx,egy,egz: components of photon momentum
     vz: target vertex of interaction
     q2: true invariant mass of photon
     phir: overall phi rotation angle in local sector system
     
## ntuple 11
       F:         0.200000E+00 
       D:         0.500000E-02 
       G:         0.300000E+01 
       TL:        0.400000E+01
       TR:        0.600000E+00 
       VX:        0.000000E+00
       VY:        0.000000E+00
       VZ:        0.000000E+00
       BE:        0.106000E+02
       EMN:       0.250000E+01
       SMN:       0.500000E+01
       SMX:       0.200000E+02
       P:         0.000000E+00
       M:         0.100000E+06
       C:         0.200000E-04
       OPST0:     1.000000E+00
       OPST1:     1.000000E+00
       OPTFI:     0.000000E+00
       IDMOF:     0.123000E+03

     
## ntuple 12  
     TELAP: Elapsed CPU time
     NEVEN: Number of generated events
     MCALL: mcall_max
     SIGI: direct integrated cross section
     SIGM: monte-carlo integrated cross section
     BMTIM: Total beam time in seconds at lum=10**34
     
## output summary files
* [elast_gen.out](https://github.com/forcar/elastgen/blob/master/elast_gen.out)
* [elast_gen.sum](https://github.com/forcar/elastgen/blob/master/elast_gen.sum)  
