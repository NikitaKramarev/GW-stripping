# GW-stripping
Function "stripping" in stripping.py calculates final stages of evolution of NS-NS binaries

## input variables in "stripping":
m1 - mass of neutron star (NS) in M_solar

m2 - mass of low-mass NS in M_solar

r1 - radius of NS in km (const in this approximation)

a0  - initial distance between NSs in km

## additional parameters:
rel = 'on' or 'off' - relativistic correction for the Roche lobe (see Ratcovic et al.)

eta = (K0 * L^2)^(1/3) - auxiliary parameter in MeV (see Sotani et al.)

max_time of calculations in s
    
## output files:
 
'stripping_dist_mass.dat' - [time in s; distance between NSs in km; q=m2/M; m2 in M_solar]
 
'stripping_rad.dat' - [time is s; GW luminosity in erg/s; nutrino luminosity in erq/s]
 
## Additional functions in stripping.py: "mass" and "roche".
 
## See description_rus.pdf and articles for details.
