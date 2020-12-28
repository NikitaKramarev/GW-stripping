import gw_stripping
import math
import fire
import sys
import numpy as np
import pandas


def er(msg):
    raise Exception(f"Fatal error:{msg}")

def write_res(res, out):
    df = pandas.DataFrame(data=res)
    s = df.to_csv(header=None, index=False)
    out.write(s)



def run(m1, m2, r1, a0, eta, rel, max_time):
    M = m1+m2                             # total mass
    q = m2/M                              # mass ratio (not according to Bisicalo!)
    x1 = 2.953286590373                   # 2*G*M_solar/km*c**2
    t0 = math.sqrt(2*a0/(x1*M))
    t0 = t0*a0/(299792.458)               # the characteristic time scale in s
    E0 = x1*M**2/(2*a0)*1.7871746729e+54  # the characteristic energy in erg
    L_gw = (2/5)*(x1/a0)**4
    L_gw = (L_gw/a0)*5.357734049e+59      # the characteristic GW luminosity in erg/s
    alpha_G = (2*x1*M/a0)**(5/2)/5        # parameter
    eps = 1e-8

    f, f_q, f_a = gw_stripping.roche(q, rel, M, a0)
    r2, m2_r2, r2_m2 = gw_stripping.mass_to_radius(m2, eta)

    if r2 > a0*f:
        er('error: start R2 > R_roche, increase start distance a0!')

    if (rel == 'off'):
        # start (dimensionless) distance of mass transfer
        a1 = r2/(a0*f)
        t1 = (1-a1**4)/(8*q*(1-q)*alpha_G)
    else:
        x1 = 10
        x2 = 1
        while abs((x2-x1)/x1) > eps:
            x1 = x2
            a = x1*a0
            f, f_q, f_a = gw_stripping.roche(q, rel, M, a)
            x2 = 1-a*f/r2
            x2 = x2*r2/(a0*(f+f_a))+x1
        a1 = x2  # start (dimensionless) distance of mass transfer (rel = 'on')
        t1 = (1-a1**4)/(8*q*(1-q)*alpha_G)

    N1 = 1000
    t = 0
    x2 = 1
    LG = M**5*q**2*(1-q)**2/x2**5
    LG = LG*L_gw

    tau1 = t1/N1
    res_init = [-t1*t0, a0, q, m2, LG]

    t_array = np.linspace(0, t1, N1)

    x2_array = np.empty(N1)
    for i, x in enumerate(t_array): x2_array[i] = ((1-(8*q*(1-q)*alpha_G*x))**(1/4)) * a0
    x1 = x2_array[-1] / a0

    time_array = np.empty(N1)
    for i, x in enumerate(t_array): time_array[i] = (x-t1)*t0
    
    q_array = np.full(N1, q)
    m2_array = np.full(N1, m2)

    LG_array = np.empty(N1)
    for i, x in enumerate(x2_array): LG_array[i] = (M**5*q**2*(1-q)**2/(x/a0)**5) * L_gw 

    res1 = list(zip(time_array, x2_array, q_array, m2_array, LG_array))

    time2_array = np.empty(N1 + 1)
    for i, x in enumerate(t_array): time2_array[i+1] = (x/2 - t1) * t0

    L_nu = np.empty(N1)

    res2 = list(zip(time2_array, L_nu[1:-1]))

    ###########################################################################
    # Stage 2: R2=R_roche, stable mass transfer

    t = t1
    tau2 = tau1
    chi = 0.55                            # parameter of the scheme

    a = x1*a0
    f, f_q, f_a = gw_stripping.roche(q, rel, M, a)
    r2 = a0*x1*f
    m2, m2_r2, r2_m2 = gw_stripping.radius_to_mass(r2, eta)
    q1 = q
    stab_corr = 1+f_a/f
    stab = f_q*q1/f-2*(1-2*q1)*stab_corr/(1-q1)
    corr = 0.005

    while r2_m2 > (stab+corr):              # stability testing

        delta = 10
        q2 = q1*(1-1e-6)
        x2 = x1*(1+1e-6)

        # solving the system of nonlinear equations f1 and f2
        while delta > eps:
            a = x2*a0
            f, f_q, f_a = gw_stripping.roche(q2, rel, M, a)
            r2 = a0*x2*f
            m2, m2_r2, r2_m2 = gw_stripping.radius_to_mass(r2, eta)

            q_2 = 1-q2
            f1 = 1/(q2*q_2*math.sqrt(x2))
            f1 = f1-1/(q1*(1-q1)*math.sqrt(x1))
            f1 = f1-alpha_G*tau2*(chi/x2**(9/2)+(1-chi)/x1**(9/2))
            f2 = q2-m2/M

            f1_x2 = -1/(2*q2*q_2*x2**(3/2))
            f1_x2 = f1_x2+9*alpha_G*tau2*chi/(2*x2**(11/2))
            f1_q2 = (2*q2-1)/(math.sqrt(x2)*(q2*q_2)**2)
            f2_x2 = -m2_r2*m2/(M*x2)
            f2_x2 = f2_x2-m2_r2*m2*f_a/(M*x2*f)
            f2_q2 = 1-m2_r2*m2*f_q/(M*f)

            delta_x = f2*f1_q2-f1*f2_q2
            delta_x = delta_x/(f1_x2*f2_q2-f2_x2*f1_q2)
            delta_q = (-f2-f2_x2*delta_x)/f2_q2
            if abs(delta_x/x2) > 0.1 or abs(delta_q/q2) > 0.1:
                k = min(abs(delta_x/x2), abs(delta_q/q2), 0.1)
            else:
                k = 1
            x2 = x2+delta_x*k
            q2 = q2+delta_q*k

            delta = math.sqrt(delta_x**2+delta_q**2)
            delta = delta/math.sqrt(x2**2+q2**2)

        t = t+tau2

        # Calculations of luminosities
        LG = M**5*q2**2*(1-q2)**2/x2**5
        LG = LG*L_gw
        L_nu = (1-(q1+q2)/2)*(q1-q2)/tau2
        L_nu = L_nu*(E0*a0/(t0*r1))

        line = ((t-t1)*t0, x2*a0, q2, q2*M, LG)
        res1.append(list(line))

        line = ((t-tau2/2-t1)*t0, L_nu)
        res2.append(list(line))

        if (abs((q2-q1)/q1) < 1.e-4 and abs((x2-x1)/x1) < 1.e-4):
            tau2 = tau2*2
        if (abs((q2-q1)/q1) > 5.e-3 or abs((x2-x1)/x1) > 5.e-3):
            tau2 = tau2/2

        x1 = x2
        q1 = q2

        a = x2*a0
        f, f_q, f_a = gw_stripping.roche(q2, rel, M, a)
        r2 = a0*x2*f
        m2, m2_r2, r2_m2 = gw_stripping.radius_to_mass(r2, eta)
        stab_corr = 1+f_a/f
        stab = f_q*q2/f-2*(1-2*q2)*stab_corr/(1-q2)

        if ((t-t1)*t0 > max_time):
            print(
                'Calculations were stopped because time of stable mass transfer > ', max_time, 'sec !', file=sys.stderr)
            break

    return res1, res2


out1_path_default = "stripping_dist_mass.dat"
out2_path_default = "stripping_rad.dat"
def main(m1=1.4, m2=0.3, r1=10, a0=100, eta=110, rel="off",
 max_time=1, out1_path=out1_path_default, out2_path=out2_path_default):
    """
    Main function
    Final evolution of close neutron star binary, stripping model

    input variables:
    m1(float): mass of neutron star (NS) in M_solar
    m2(float): mass of low-mass NS in M_solar
    r1(float): radius of NS in km
    a0(float): initial distance between NSs in km
    out1_path: stripping dist mass output file
    out2_path: stripping rad output file

    additional parameters:
    eta(float): (K0*L**2)**(1/3), auxiliary parameter in MeV (see Sotani et al.)
    rel(str): 'on' or 'off', relativistic correction for the Roche lobe
    max_time(float): stopping time in s

    output files:
    'stripping.dat' - [time in s; distance between NSs in km; q=m2/M; m2 in M_solar; GW luminosity in erg/s]
    'stripping_rad.dat' - [time in s; nutrino luminosity in erq/s]

    note:
    see file description_rus.pdf for details
    """

    if not (rel in ("on", "off")):
        er("wrong rel")

    if m1 < m2:
        er('m1 must be more than m2!')

    if eta > 200:
        er('error: eta must be less than 200!')

    if eta < 60:
        er("error: eta must be more than 60!")

    if not (m2 > 0.1 and m2 < 0.8):
        er('error: m2 must be in (0.1, 0.8)!')

    if m2 < 0: 
        er('error: m2 must be positive!')

    if m2 > 1: 
        er('warning, mass of m1 is more than 1 M_sol so approximation of Sotani may be not valid!')
    
    if m1 > 2.3:
        print('Warning, m1 have a mass of BH!')

    
    with open(out1_path, 'w') as out1:
        with open(out2_path, 'w') as out2:
            res1, res2 = run(m1=m1, m2=m2, r1=r1, a0=a0, eta=eta, rel=rel, max_time=max_time)
            write_res(res1, out1)
            write_res(res2, out2)



if __name__ == "__main__":
    # m1=1.4
    # m2=0.3
    # r1=10
    # a0=100
    # eta=110
    # rel="on"
    # max_time=1
    # run(m1=m1, m2=m2, r1=r1, a0=a0, eta=eta, rel=rel, max_time=max_time)
    fire.Fire(main)