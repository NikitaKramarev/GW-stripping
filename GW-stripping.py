import math
#import numpy as np
#import matplotlib.pyplot as plt
############################################################################################
# function 1
# Approximations to the Radii of Roche Lobes
#
# input variables:
# input variables:
# q = m2/M
# M = m1 + m2 - total mass in M_solar
# a - distance between NSs in km
# rel = 'on' or 'off' - relativistic correction 
#
# output variables:
# f(q) = R_roche/a
# f_q = d(f(q))/dq
# f_a = a*d(f(a))/da (if rel = 'on')

def roche(q, rel, M, a):
    qs = q/(1-q)                       # qs = m2/m1
    aq = 0.49
    bq = 0.6
    x1 = qs**(2/3)
    x2 = 1+qs**(1/3)
    if rel == 'off':                   # nonrelativistic approximation (according to Eggleton, 1982)
        f = (aq*x1)
        f = f/(bq*x1+math.log(x2))

        f_q = (2*aq*qs**(-1/3)*math.log(x2)/3-aq/(3*x2))
        f_q = f_q/(bq*x1+math.log(x2))**2
        f_q = f_q/(1-q)**2
        f_a = 0

    elif rel == 'on':                  # correction that stems from post-Newtonian effects (according to Ratcovic et al.)
        z = 2.953286590373*M/a
        ac = 1.951
        bc = 1.812
        
        qrel = (aq*x1)
        qrel = qrel/(bq*x1+math.log(x2))
        crel = 1+z*(ac*qs**(1/5)-bc)
        qrel_q = (2*aq*qs**(-1/3)*math.log(x2)/3-aq/(3*x2))
        qrel_q = qrel_q/(bq*x1+math.log(x2))**2
        crel_q = z*ac*qs**(-4/5)/5

        f = qrel*crel
        f_q = qrel*crel_q+qrel_q*crel
        f_q = f_q/(1-q)**2
        f_a = -qrel*(ac*qs**(1/5)-bc)*z

    else:
        print('Error: character "rel" must be "on" or "off"!')
        raise SystemExit(1)
        
    return f, f_q, f_a

############################################################################################
# function 2
# Mass-radius formulas for low-mass neutron star in NS-NS or BH-NS coalescing binaries
# according to Sotani et al.
#
# input variables:
# r2 - radius in km or
# m2 - mass in M_solar
# m2_or_r2 = 'm2' or 'r2'
# eta = (K0*L**2)**(1/3) - auxiliary parameter in MeV (see Sotani et al.)
#
# output variables:
# m2 or r2
# m2_r2 = d(ln(m2))/d(ln(r2))
# r2_m2 = d(ln(r2))/d(ln(m2))

def mass(input_var, m2_or_r2, eta):
    eta = eta/100                      # normalisation 
    x1 = 2.953286590373                # 2*G*M_solar/km*c**2
    a1 = 0.279-0.235*eta
    b1 = -0.82+1.25*eta
    a2 = 0.0255-0.012*eta
    b2 = 0.108*eta-0.0619

    if m2_or_r2 == 'r2':
        
        r2 = input_var
        m2 = 0.2                       # initial approach for mass in M_solar
        delta = 10
        while abs(delta) > 1e-11:      # search of m2
            c1 = 0.371-0.593*eta-m2
            uc = -b1+math.sqrt(b1**2-4*a1*c1)
            uc = uc/(2*a1)             # uc = rho_c/rho_nucl
            c2 = 0.00859-0.0429*eta
            z = a2*uc**2+b2*uc+c2      # the gravitational redshift
            z1 = z+1
            z2 = z+2
            r2_2 = (z+1)**2/(z**2+2*z)
            r2_2 = r2_2*x1*m2          # in km (!!!)

            r2_m2 = x1*m2*2*z1*(z*z2-z1**2)/(z**2*z2**2)
            r2_m2 = r2_m2*(2*a2*uc+b2)
            c1 = 0.371-0.593*eta-m2
            r2_m2 = r2_m2/math.sqrt(b1**2-4*a1*c1)
            r2_m2 = r2_m2 + x1*z1**2/(z*z2)

            delta = (r2_2-r2)/(r2_m2*m2)
            k=1
            if (abs(delta) > 0.1): k=0.1
            m2 = m2*(1-delta*k)         # m2 in M_solar 
            
        c1 = 0.371-0.593*eta-m2
        uc = -b1+math.sqrt(b1**2-4*a1*c1)
        uc = uc/(2*a1)
        c2 = 0.00859-0.0429*eta
        z = a2*uc**2+b2*uc+c2
        z1 = z+1
        z2 = z+2
        output_var = m2
        
    elif m2_or_r2 == 'm2':
        
        m2 = input_var
        c1 = 0.371-0.593*eta-m2
        uc = -b1+math.sqrt(b1**2-4*a1*c1)
        uc = uc/(2*a1)
        c2 = 0.00859-0.0429*eta
        z = a2*uc**2+b2*uc+c2
        z1 = z+1
        z2 = z+2
        r2 = x1*m2*z1**2/(z*z2)        # r2 in km 
        output_var = r2

    else:
        print('Error: character "m2_or_r2" must be "r2" or "m2"!')
        raise SystemExit(1) 
            
    r2_m2 = x1*m2*2*z1*(z*z2-z1**2)/(z**2*z2**2)
    r2_m2 = r2_m2*(2*a2*uc+b2)
    c1 = 0.371-0.593*eta-m2
    r2_m2 = r2_m2/math.sqrt(b1**2-4*a1*c1)
    r2_m2 = r2_m2 + x1*z1**2/(z*z2)
    r2_m2 = r2_m2*m2/r2                # d(ln(r2))/d(ln(m2))
    m2_r2 = 1/r2_m2                    # d(ln(m2))/d(ln(r2))    

    if uc > 2:
        print('Warning, uc =', uc, ' (rho_c > 2*rho_nucl), so mass-radius approximation of Sotani may be not valid!')
        input('Press to continue...')
    elif uc < 0.9:
        print('Warning, uc =', uc, ' (rho_c < 0.9*rho_nucl), so mass-radius approximation of Sotani may be not valid!')
        input('Press to continue...')
            
    return output_var, m2_r2, r2_m2

############################################################################################
# Main function
# Stage 1: R2<R_roche
#
# input variables:
# m1 - mass of neutron star (NS) in M_solar
# m2 - mass of low-mass NS in M_solar
# r1 - radius of NS in km
# a0  - initial distance between NSs in km
#
# additional parameters:
# eta = (K0*L**2)**(1/3) - auxiliary parameter in MeV (see Sotani et al.)
# rel = 'on' or 'off' - relativistic correction for the Roche lobe
# max_time in s
#
# output files:
# 'stripping_dist_mass.dat' - [time in s; distance between NSs in km; q=m2/M; m2 in M_solar]
# 'stripping_rad.dat' - [time is s; GW luminosity in erg/s; nutrino luminosity in erq/s]

m1  = 1.4 
m2  = 0.3 
r1  = 10
a0  = 100
eta = 110 
rel = 'on'
max_time = 3

M  = m1+m2                            # total mass
q  = m2/M                             # mass ratio (not according to Bisicalo!)
x1 = 2.953286590373                   # 2*G*M_solar/km*c**2
t0 = math.sqrt(2*a0/(x1*M))           
t0 = t0*a0/(299792.458)               # the characteristic time scale in s
E0 = x1*M**2/(2*a0)*1.7871746729e+54  # the characteristic energy in erg
alpha_G = (2*x1*M/a0)**(5/2)/5        # parameter
eps = 1e-8

f, f_q, f_a = roche(q, rel, M, a0)
r2, m2_r2, r2_m2 = mass(m2, 'm2', eta)

if m1 < m2:
    print('Error: m1 must be more than m2!')
    raise SystemExit(1)
if r2 > a0*f:
    print('Error: start R2 > R_roche, increase start distance a0!')
    raise SystemExit(1)

if (rel == 'off'):
    a1 = r2/(a0*f)                    # start (dimensionless) distance of mass transfer
    t1 = (1-a1**4)/(8*q*(1-q)*alpha_G) 
else:
    x1 = 10
    x2 = 1
    while abs((x2-x1)/x1)>eps:
        x1 = x2
        a = x1*a0
        f, f_q, f_a = roche(q, rel, M, a)
        x2 = 1-a*f/r2
        x2 = x2*r2/(a0*(f+f_a))+x1
    a1 = x2                           #start (dimensionless) distance of mass transfer (rel = 'on')
    t1 = (1-a1**4)/(8*q*(1-q)*alpha_G) 

N1 = 1000
t = 0

file1 = open('stripping_dist_mass.dat', 'w')
file2 = open('stripping_rad.dat', 'w')
line = (-t1*t0,a0,q,m2)
for s in line:
    file1.write('%20.8e' % s)         # Writing part
file1.write('\n')     

tau1 = t1/N1
q_1 = 1-q
x1 = 1
L_nu = 0
for i in range(N1):                   # Writing part and calculations of LG luminosity 
    t  = t+tau1
    x2  = (1-(8*q*(1-q)*alpha_G*t))**(1/4)
    LG = q*q_1*(1/x2-1/x1)/(2*tau1)
    LG = LG*(E0/t0)
    line = ((t-t1)*t0,x2*a0,q,m2)
    for s in line:
        file1.write('%20.8e' % s)       
    file1.write('\n')    
    line = ((t-tau1/2-t1)*t0,LG,L_nu)
    for s in line:
        file2.write('%20.8e' % s)       
    file2.write('\n') 
    x1 = x2   
############################################################################################
# Stage 2: R2=R_roche, stable mass transfer

tau2 = tau1
chi   = 0.55                            # parameter of the scheme

a = x1*a0
f, f_q, f_a = roche(q, rel, M, a)
r2 = a0*x1*f
m2, m2_r2, r2_m2 = mass(r2, 'r2', eta)
q1 = q
stab = f_q*q1/f-2*(1-2*q1)/(1-q1)
corr = 0.005

while r2_m2 > (stab+corr):              # stability testing

    delta = 10
    q2 = q1*(1-1e-6)
    x2 = x1*(1+1e-6)
    
    while delta > eps:                  # solving the system of nonlinear equations f1 and f2        

        a = x2*a0       
        f, f_q, f_a = roche(q2, rel, M, a)
        r2 = a0*x2*f
        m2, m2_r2, r2_m2 = mass(r2, 'r2', eta)

        q_2 = 1-q2
        f1  = 1/(q2*q_2*math.sqrt(x2))
        f1  = f1-1/(q1*(1-q1)*math.sqrt(x1))
        f1  = f1-alpha_G*tau2*(chi/x2**(9/2)+(1-chi)/x1**(9/2))
        f2  = q2-m2/M

        f1_x2 = -1/(2*q2*q_2*x2**(3/2))
        f1_x2 = f1_x2+9*alpha_G*tau2*chi/(2*x2**(11/2))
        f1_q2 = (2*q2-1)/(math.sqrt(x2)*(q2*q_2)**2)
        f2_x2 = -m2_r2*m2/(M*x2)
        f2_x2 = f2_x2-m2_r2*m2*f_a/(M*x2*f)  
        f2_q2 = 1-m2_r2*m2*f_q/(M*f)

        delta_x = f2*f1_q2-f1*f2_q2
        delta_x = delta_x/(f1_x2*f2_q2-f2_x2*f1_q2)
        delta_q = (-f2-f2_x2*delta_x)/f2_q2
        if abs(delta_x/x2)>0.1 or abs(delta_q/q2)>0.1:
            k = min(abs(delta_x/x2),abs(delta_q/q2),0.1)
        else:
            k=1
        x2 = x2+delta_x*k
        q2 = q2+delta_q*k

        delta   = math.sqrt(delta_x**2+delta_q**2)
        delta   = delta/math.sqrt(x2**2+q2**2)

    t = t+tau2

    LG = q1*(1-q1)/x1                   # Writing part and calculations of luminosities
    LG = LG-q2*(1-q2)/x2
    LG = LG/(2*tau2)*(E0/t0)
    L_nu = (1-(q1+q2)/2)*(q1-q2)/tau2
    L_nu = L_nu*(E0*a0/(t0*r1))
    line = ((t-t1)*t0,x2*a0,q2,q2*M)
    for s in line:
        file1.write('%20.8e' % s)       
    file1.write('\n') 
    line = ((t-tau2/2-t1)*t0,LG,L_nu)
    for s in line:
        file2.write('%20.8e' % s)       
    file2.write('\n')

    if (abs((q2-q1)/q1)<1.e-4 and abs((x2-x1)/x1)<1.e-4): tau2=tau2*2
    if (abs((q2-q1)/q1)>5.e-3 or abs((x2-x1)/x1)>5.e-3): tau2=tau2/2
    
    x1 = x2
    q1 = q2

    a = x2*a0
    f, f_q, f_a = roche(q2, rel, M, a)
    r2 = a0*x2*f
    m2, m2_r2, r2_m2 = mass(r2, 'r2', eta)
    stab = f_q*q2/f-2*(1-2*q2)/(1-q2)

    if ((t-t1)*t0>max_time): 
        print('Calculations were stopped because time of stable mass transfer > ', max_time,'sec !')
        break

file1.close()
file2.close()
print('Done!')

############################################################################################

    




        


    
    
