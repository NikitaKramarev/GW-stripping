import math
import sys
from scipy import optimize

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

class mass_tech:
    def mass_to_radius(self, input_var, suspendWarning=False):
        m2 = input_var
        c1 = 0.371-0.593*self.eta-m2
        uc = -self.b1+math.sqrt(self.b1**2-4*self.a1*c1)
        uc = uc/(2*self.a1)
        c2 = 0.00859-0.0429*self.eta
        z = self.a2*uc**2+self.b2*uc+c2
        z1 = z+1
        z2 = z+2
        r2 = self.x1*m2*z1**2/(z*z2)        # r2 in km
        output_var = r2

        r2_m2 = self.x1*m2*2*z1*(z*z2-z1**2)/(z**2*z2**2)
        r2_m2 = r2_m2*(2*self.a2*uc+self.b2)
        c1 = 0.371-0.593*self.eta-m2
        r2_m2 = r2_m2/math.sqrt(self.b1**2-4*self.a1*c1)
        r2_m2 = r2_m2 + self.x1*z1**2/(z*z2)
        r2_m2 = r2_m2*m2/r2                # d(ln(r2))/d(ln(m2))
        m2_r2 = 1/r2_m2                    # d(ln(m2))/d(ln(r2))

        if (suspendWarning == False):
            self.check(uc)

        return output_var, m2_r2, r2_m2


    def radius_to_mass(self, input_var):
        r2 = input_var

        sol = optimize.root_scalar(lambda x: input_var - list(self.mass_to_radius(x, True))[0], bracket=[0.1, 0.8], method='brentq')
        m2 = sol.root

        c1 = 0.371-0.593*self.eta-m2
        uc = -self.b1+math.sqrt(self.b1**2-4*self.a1*c1)
        uc = uc/(2*self.a1)
        c2 = 0.00859-0.0429*self.eta
        z = self.a2*uc**2+self.b2*uc+c2
        z1 = z+1
        z2 = z+2
        output_var = m2

        r2_m2 = self.x1*m2*2*z1*(z*z2-z1**2)/(z**2*z2**2)
        r2_m2 = r2_m2*(2*self.a2*uc+self.b2)
        c1 = 0.371-0.593*self.eta-m2
        r2_m2 = r2_m2/math.sqrt(self.b1**2-4*self.a1*c1)
        r2_m2 = r2_m2 + self.x1*z1**2/(z*z2)
        r2_m2 = r2_m2*m2/r2                # d(ln(r2))/d(ln(m2))
        m2_r2 = 1/r2_m2                    # d(ln(m2))/d(ln(r2))

        self.check(uc)

        return output_var, m2_r2, r2_m2


    def check(self, uc):
        if uc > 2:
            print('Warning, uc =', uc,
                ' (rho_c > 2*rho_nucl), so mass-radius approximation of Sotani may be not valid!', file=sys.stderr)
        elif uc < 0.75:
            print('Warning, uc =', uc,
                ' (rho_c < 0.75*rho_nucl), so mass-radius approximation of Sotani may be not valid!', file=sys.stderr)

    def __init__(self, eta):
        self.eta = eta/100                      # normalisation
        self.x1 = 2.953286590373                # 2*G*M_solar/km*c**2
        self.a1 = 0.279-0.235*self.eta
        self.b1 = -0.82+1.25*self.eta
        self.a2 = 0.0255-0.012*self.eta
        self.b2 = 0.108*self.eta-0.0619

def mass_to_radius(inp, eta):
    """ Mass-radius formulas for low-mass neutron star in NS-NS or BH-NS coalescing binaries
    according to Sotani et al. Mass to radius.

    input variables:
        input_var(float): radius in km (r2) or mass in M_solar (m2)
        eta(float): (K0*L**2)**(1/3), auxiliary parameter in MeV (see Sotani et al.)

    output variables:
        output_var(float): mass in M_solar (m2) or radius in km (r2)
        m2_r2(float): d(ln(m2))/d(ln(r2)), logarithmic derivative
        r2_m2(float): d(ln(r2))/d(ln(m2)), logarithmic derivative
    """
    m = mass_tech(eta)
    return m.mass_to_radius(inp)

def radius_to_mass(inp, eta):
    """ Mass-radius formulas for low-mass neutron star in NS-NS or BH-NS coalescing binaries
    according to Sotani et al. Radius to mass

    input variables:
        input_var(float): radius in km (r2)
        eta(float): (K0*L**2)**(1/3), auxiliary parameter in MeV (see Sotani et al.)

    output variables:
        output_var(float): mass in M_solar (m2) or radius in km (r2)
        m2_r2(float): d(ln(m2))/d(ln(r2)), logarithmic derivative
        r2_m2(float): d(ln(r2))/d(ln(m2)), logarithmic derivative
    """
    m = mass_tech(eta)
    return m.radius_to_mass(inp)
