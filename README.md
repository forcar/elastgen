# elastgen
c     Original author: Ralph Minehart (University of Virginia)
c     Modifications by K. Joo and L.C. Smith (UVA)
c     See http://clasweb.jlab.org/cgi-bin/cvsweb/cvsweb.cgi/packages/aao/elast_gen

c     This program makes ntuples and BOS event files for (e,e'p) elastic scattering
c     with the elastic peak modified by vertex corrections according to the Mo-Tsai 
c     formula (II.6) and internal bremmstrahlung (radiative tails) calculated from 
c     a Monte-Carlo sampling of formula (B.4) and (B.5) from RMP, Vol. 41, 205 (1969). 

c     Ntuple 10 contains:
c     ES: incident e- energy 
c     EP: scattered e- energy
c     THETE: scattered e- theta
c     W:  true hadronic invariant mass
c     PPRX, PPRY, PPRZ: components of p+ momentum
c     EPROT: proton energy  
c     THETP: proton theta
c     PHIPR: proton phi
c     MM: missing mass
c     CSTHK,PHIK: photon angles relative to the q vector
c     EG: photon energy
c     EGX,EGY,EGZ: components of photon momentum
c     VX,VY,VZ: target vertex of interaction
c     Q2: true invariant mass of photon

c     Ntuple 11 contains:
c     Input arguments after file4

c     Ntuple 12 contains:
c     TELAP: Elapsed CPU time
c     NEVEN: Number of generated events
c     MCALL: mcall_max
c     SIGI: direct integrated cross section
c     SIGM: monte-carlo integrated cross section
c     BMTIM: Total beam time in seconds at lum=10**34
      
c     Some tricks are employed to make the monte-carlo concentrate
c     on the regions where the cross section is large.  Also the
c     code attempts to find the maximum value of the radiative integrand
c     to normalize the monte-carlo sampling (sigr_max).  Note that 
c     sigr_max may be underestimated when incident e- straggling is turned
c     on.  The code attempts to compensate for this by oversampling those
c     events with large incident e- energy loss (mcall>1), which may cause 
c     spurious statistical spikes when f=1 and delta=0. For that choice of
c     options, input e- straggling should be commented out.  

c     This version of the program makes only one new file and then quits.
