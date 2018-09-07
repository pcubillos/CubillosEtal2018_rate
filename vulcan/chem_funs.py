import os
import numpy as np

root = os.path.realpath(os.path.dirname(__file__))

#spec_list = ['H', 'H2O', 'OH', 'H2', 'O', 'CH', 'C', 'CH2', 'CH3', 'CH4', 'C2', 'C2H2', 'C2H', 'C2H3', 'C2H4', 'C2H5', 'C2H6', 'CO', 'CO2', 'CH2OH', 'H2CO', 'HCO', 'CH3O', 'CH3OH', 'CH3CO', 'O2', 'H2CCO', 'HCCO', 'He']
spec_list = ['H2O', 'H2', 'CH4', 'C2H2', 'C2H4', 'CO', 'CO2']


# the data of 'H2CO' is from Brucat's 2015
#C2H NASA 9 new from Brucat
nasa9 = {}
nasa9['He','low'], nasa9['He','high'] = [0]*10, [0]*10
for i in [ _ for _ in spec_list + ['N2', 'NH3', 'HCN'] if _ != 'He']:
    nasa9[i] = np.loadtxt(root+'/thermo_NASA9/' + str(i) + '.txt')
    nasa9[i] = nasa9[i].flatten()
    nasa9[i,'low'] = nasa9[i][0:10]
    nasa9[i,'high'] = nasa9[i][10:20]

#H/RT
def h_RT(T,a):
    return -a[0]/T**2 + a[1]*np.log(T)/T + a[2] + a[3]*T/2. + a[4]*T**2/3. + a[5]*T**3/4. + a[6]*T**4/5. + a[8]/T

#s/R
def s_R(T,a):
    return -a[0]/T**2/2. -a[1]/T + a[2]*np.log(T) + a[3]*T + a[4]*T**2/2. + a[5]*T**3/3. + a[6]*T**4/4. + a[9]

#g/RT (non-dimensional)
def g_RT(T,a_low,a_high):
    #print a_low[:2]
    if np.any(np.logical_or(T < 200, T > 6000)):
        print 'T exceeds the valid range.'
        gi = 'T out of data range!'

    gi = (T < 1000)*(h_RT(T, a_low)-s_R(T, a_low)) + (T >= 1000)*(h_RT(T, a_high)-s_R(T, a_high))

    return gi

def gibbs_sp(i,T):
    gi = g_RT(T,nasa9[i,'low'],nasa9[i,'high'])
    return gi
