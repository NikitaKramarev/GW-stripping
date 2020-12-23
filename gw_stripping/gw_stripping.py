import math


def roche_nonrelativistic(q):
    """
    Approximations to the Radii of Roche Lobes according to Eggleton or Ratcovic

    input variables:  
    q(float): m2/M, mass ratio  

    output variables:  
    f(float): R_roche/a  
    f_q(float): d(f(q))/dq, derivative
    """

    if (q == 0):
        return 0, 0

    qs = q/(1-q)                       # qs = m2/m1
    aq = 0.49
    bq = 0.6
    x1 = qs**(2/3)
    x2 = 1+qs**(1/3)

    # nonrelativistic approximation (according to Eggleton, 1982)
    f = (aq*x1)
    f = f/(bq*x1+math.log(x2))

    f_q = (2*aq*qs**(-1/3)*math.log(x2)/3-aq/(3*x2))
    f_q = f_q/(bq*x1+math.log(x2))**2
    f_q = f_q/(1-q)**2

    return f, f_q


def roche_relativistic(q, M, a):
    """ 
    Approximations to the Radii of Roche Lobes according to Eggleton or Ratcovic

    input variables:  
    q(float): m2/M, mass ratio  
    M(float): m1 + m2, total mass of the binary NSs system in M_solar  
    a(float): distance between NSs in km  

    output variables:
    f(float): R_roche/a  
    f_q(float): d(f(q))/dq, derivative  
    f_a(float): a*d(f(a))/da 
  """
    qs = q/(1-q)                       # qs = m2/m1
    aq = 0.49
    bq = 0.6
    x1 = qs**(2/3)
    x2 = 1+qs**(1/3)

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

    return f, f_q, f_a


def roche(q, rel, M=0, a=0):
    if (rel == "on"):
        return roche_relativistic(q, M, a)
    else:
        r1, r2, = roche_nonrelativistic(q)
        return r1, r2, 0.0


def mass(input_var, m2_or_r2, eta):
    """ 
Mass-radius formulas for low-mass neutron star in NS-NS or BH-NS coalescing binaries
according to Sotani et al.

  input variables:
    input_var(float): radius in km (r2) or mass in M_solar (m2)
    m2_or_r2(str): 'm2' or 'r2', type of variable
    eta(float): (K0*L**2)**(1/3), auxiliary parameter in MeV (see Sotani et al.)

  output variables:
    output_var(float): mass in M_solar (m2) or radius in km (r2)
    m2_r2(float): d(ln(m2))/d(ln(r2)), logarithmic derivative
    r2_m2(float): d(ln(r2))/d(ln(m2)), logarithmic derivative

  raises: SystemExit
  """
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
            r2_2 = r2_2*x1*m2          # in km

            r2_m2 = x1*m2*2*z1*(z*z2-z1**2)/(z**2*z2**2)
            r2_m2 = r2_m2*(2*a2*uc+b2)
            c1 = 0.371-0.593*eta-m2
            r2_m2 = r2_m2/math.sqrt(b1**2-4*a1*c1)
            r2_m2 = r2_m2 + x1*z1**2/(z*z2)

            delta = (r2_2-r2)/(r2_m2*m2)
            k = 1
            if (abs(delta) > 0.1):
                k = 0.1
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
        print('Warning, uc =', uc,
              ' (rho_c > 2*rho_nucl), so mass-radius approximation of Sotani may be not valid!')
        # input('Press to continue...')
    elif uc < 0.75:
        print('Warning, uc =', uc,
              ' (rho_c < 0.75*rho_nucl), so mass-radius approximation of Sotani may be not valid!')
        # input('Press to continue...')

    return output_var, m2_r2, r2_m2
