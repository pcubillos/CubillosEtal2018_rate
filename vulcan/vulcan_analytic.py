import numpy as np
from numpy import polynomial

from chem_funs import gibbs_sp

n_eq, mix_eq = {}, {}

kk1 = lambda T, pbar: np.exp( -( gibbs_sp('CO',T) + 3*gibbs_sp('H2',T)
                    -1*gibbs_sp('CH4',T)-1*gibbs_sp('H2O',T) ) ) * (pbar)**-2
kk2 = lambda T: np.exp( -( gibbs_sp('CO',T) + gibbs_sp('H2O',T)
                        -1*gibbs_sp('CO2',T) -1*gibbs_sp('H2',T) ) )
kk3 = lambda T, pbar: np.exp( -( gibbs_sp('C2H2',T) + 3*gibbs_sp('H2',T)
                              -2*gibbs_sp('CH4',T) ) ) * (pbar)**-2
kk4 = lambda T, pbar: np.exp( -( gibbs_sp('C2H2',T) + gibbs_sp('H2',T)
                                -gibbs_sp('C2H4',T) ) ) * (pbar)**-1
kk5 = lambda T, pbar: np.exp( -( gibbs_sp('N2',T) + 3*gibbs_sp('H2',T)
                              -2*gibbs_sp('NH3',T) ) ) * (pbar)**-2
kk6 = lambda T, pbar: np.exp( -( gibbs_sp('HCN',T) + 3*gibbs_sp('H2',T)
                        -gibbs_sp('NH3',T) -gibbs_sp('CH4',T) ) ) * (pbar)**-2

sp_ini = ["H2O", "CO", "CO2", "CH4", "C2H2", "C2H4", "HCN", "NH3", "N2"]
He_H = 0.09


def run_analytic_vulcan(temp, press, n_o, n_c, n_n):
  """
  Extract from VULCAN code: https://github.com/exoclime/VULCAN
  """
  nlayers = len(temp)
  q = np.zeros((nlayers, len(sp_ini)+2))
  for i in np.arange(nlayers):
    tp, pp = temp[i], press[i]
    k1 = kk1(tp,pp)
    k5 = kk5(tp,pp)
    k6 = kk6(tp,pp)
    d2 = 1.0 + 2.0*k1*(n_o+n_c)
    a0 = 256.0*(k1**3)*k5*(n_o**3)*n_c*n_c
    a1 = 32.0*((k1*n_o)**2)*n_c*( k6 - 4.0*k5*( d2 + k1*n_c ) )
    a2 = 16.0*k1*k5*n_o*( 8.0*k1*k1*n_o*n_c + d2*d2 + 4.0*k1*d2*n_c ) + 8.0*k1*k6*n_o*( 2.0*k6*(n_c-n_n) - 2.0*k1*n_c - d2 )
    a3 = -8.0*k1*k5*( 4.0*k1*d2*n_o + 8.0*k1*k1*n_o*n_c + d2*d2 ) + 4.0*k1*k6*( 2.0*k1*n_o + d2 ) + 4.0*k6*k6*( 2.0*k1*n_n - d2 )
    a4 = 16.0*k1*k1*k5*( k1*n_o + d2 ) + 4.0*k1*k6*(k6-k1)
    a5 = -8.0*(k1**3)*k5
    result = polynomial.polynomial.polyroots([a0,a1,a2,a3,a4,a5])
    result = result[result.real > 0.0]
    result = result[result.real < 2.0*n_o]
    try: # multiple positive realroots
        result = result[0]#.real
    except: # for only single positive real root
        pass
    n_eq['CO'] = result.real
    # for no root found, meaning the mixing ratio of CO is extreamly low
    if not n_eq['CO']:
        n_eq['CO']= 0.
        print ('No equilibrium solution found for CO! Setting to zero.')
    k2 = kk2(tp)
    k3 = kk3(tp,pp)
    k4 = kk4(tp,pp)
    c2 = 1.0/kk2(tp)
    n_eq['H2O'] = (2.0*n_o - n_eq['CO'])/(1.0 + 2.0*c2*n_eq['CO'])
    n_eq['CH4'] = n_eq['CO']/k1/n_eq['H2O']
    n_eq['CO2'] = n_eq['CO']*n_eq['H2O']/k2
    n_eq['C2H2'] = k3*n_eq['CH4']**2
    n_eq['C2H4'] = n_eq['C2H2']/k4
    term1 = 1. + k6*n_eq['CH4']
    n_eq['NH3'] = ( np.sqrt( term1**2 + 16.*n_n*k5 ) - term1 )/4./k5
    n_eq['HCN'] =  k6*n_eq['NH3']*n_eq['CH4']
    n_eq['N2'] = k5*n_eq['NH3']*n_eq['NH3']
    norm = sum([n_eq[sp] for sp in sp_ini]) + 1. + He_H*2. # 1. is from H2
    for sp in sp_ini: # normalizing so that the total mixing ratio equals 1
        mix_eq[sp] = n_eq[sp]/norm
    mix_eq['H2'] = 1./norm
    mix_eq['He'] = 2.*He_H/norm
    q[i] = [mix_eq["H2"],  mix_eq["He"],  mix_eq["H2O"],  mix_eq["CO"],
            mix_eq["CO2"], mix_eq["CH4"], mix_eq["C2H2"],  mix_eq["C2H4"],
            mix_eq["HCN"], mix_eq["NH3"], mix_eq["N2"]]
  return q
