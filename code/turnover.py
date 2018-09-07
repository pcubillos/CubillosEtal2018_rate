#! /usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate as si
import scipy.optimize    as so


def poly1(p, T, C, N, O):
  """
  Linear polynomial function in all variables.
  """
  mesh = np.meshgrid(C, T, N, O)
  top = p[0] + p[1]*mesh[1] + p[2]*mesh[0] + p[3]*mesh[2] + p[4]*mesh[3]
  return np.clip(top, -8, 3)


def poly4111(p, T, C, N, O):
  """
  Quartic polinomial in T, linear in CNO.
  """
  mesh = np.meshgrid(C, T, N, O)
  top = (p[0]
       + p[1]*mesh[1] + p[2]*mesh[1]**2 + p[3]*mesh[1]**3 + p[4]*mesh[1]**4
       + p[5]*mesh[0]
       + p[6]*mesh[2]
       + p[7]*mesh[3])
  return np.clip(top, -8, 3)


def poly4141(p, T, C, N, O):
  """
  Quartic polynomial in T, linear C, quartic N, linear O.
  """
  mesh = np.meshgrid(C, T, N, O)
  top = (p[ 0]
       + p[ 1]*mesh[1] + p[2]*mesh[1]**2 + p[3]*mesh[1]**3 + p[4]*mesh[1]**4
       + p[ 5]*mesh[0]
       + p[ 6]*mesh[2] + p[7]*mesh[2]**2 + p[8]*mesh[2]**3 + p[9]*mesh[2]**4
       + p[10]*mesh[3])
  return np.clip(top, -8, 3)


def poly4(p, T, C, N, O):
  """
  Quartic polynomial function in all variables.
  """
  mesh = np.meshgrid(C, T, N, O)
  top = (p[ 0]
       + p[ 1]*mesh[1] + p[ 2]*mesh[1]**2 + p[ 3]*mesh[1]**3 + p[ 4]*mesh[1]**4
       + p[ 5]*mesh[0] + p[ 6]*mesh[0]**2 + p[ 7]*mesh[0]**3 + p[ 8]*mesh[0]**4
       + p[ 9]*mesh[2] + p[10]*mesh[2]**2 + p[11]*mesh[2]**3 + p[12]*mesh[2]**4
       + p[13]*mesh[3] + p[14]*mesh[3]**2 + p[15]*mesh[3]**3 + p[16]*mesh[3]**4)
  return np.clip(top, -8, 3)


def error(p, T,C,N,O, data, model, Tmax=np.inf, Tmin=-np.inf):
  """
  Error polynomial function to optimize.
  """
  plaw = model(p,T,C,N,O)
  diff = data - plaw
  # Ignore points beyond boundaries:
  diff[data==0] = 0.0
  diff[~np.isfinite(data)] = 0.0
  # Ignore T > Tmax:
  diff[T>np.log10(Tmax)] = 0.0
  # Ignore T < Tmin:
  diff[T<np.log10(Tmin)] = 0.0
  return diff.flatten()


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read TEA data from npz file:
with np.load("TEA_grid.npz") as d:
  H2O = d["H2O"]
  CO  = d["CO"]

nlayers, nt, nz, dummy, dummy = np.shape(H2O)

# Hardcoded values (same from tea_grid.py):
press = np.logspace(-8, 3, nlayers)
temp  = np.logspace(2.305, 3.778, nt)
C = np.logspace(-6.57, -0.57, nz)
N = np.logspace(-7.17, -1.17, nz)
O = np.logspace(-6.31, -0.31, nz)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find turn-over pressure (TOP) between H2O/CO dominated atmospheres:

logp = np.log10(press)
TOP = np.zeros((nt,nz,nz,nz))

for t in np.arange(nt):
  for o in np.arange(nz):
    for c in np.arange(nz):
      if C[c] < O[o]:
        continue
      for n in np.arange(nz):
        # Ignore high-metallicity cases:
        if C[c] + O[o] + N[n] >= 0.1:
          continue
        # Index above TOP:
        if np.all(CO[:,t,c,n,o] < H2O[:,t,c,n,o]):
          TOP[t,c,n,o] = press[0]
        elif np.all(CO[:,t,c,n,o] > H2O[:,t,c,n,o]):
          TOP[t,c,n,o] = press[-1]
        else:
          # Use a spline to fine-tune interpolate TOP:
          sCO  = si.interp1d(logp, np.log10(CO [:,t,c,n,o]))
          sH2O = si.interp1d(logp, np.log10(H2O[:,t,c,n,o]))
          # 20x oversample around turn-over pressure:
          i = np.where(CO[:,t,c,n,o] > H2O[:,t,c,n,o])[0][-1]
          p = np.linspace(logp[i], logp[i+1], 20)
          # Turn-over pressure:
          i = np.where(sCO(p) > sH2O(p) )[0][-1]
          TOP[t,c,n,o] = 10**p[i]


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Fit a polynomial to TOP:

# Set all variables in log-scale:
T = temp

logT = np.log10(T)
logC = np.log10(C)
logN = np.log10(N)
logO = np.log10(O)
logTOP = np.log10(TOP)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Fit TOP with progressively increasing-order polynomials:
# (use lower-order fits as initial guess)

# Linear poly in all variables:
guess = [1.0, 0.0, 0.0, 0.0, 0.0]
fit = so.leastsq(error, guess, args=(logT,logC,logN,logO,logTOP,poly1),
                 full_output=True)
pfit1 = np.array(fit[0])

# Quartic-T + linear-CNO:
guess = np.insert(pfit1, 2, [0.0, 0.0, 0.0])
fit = so.leastsq(error, guess, args=(logT,logC,logN,logO,logTOP,poly4111),
                 full_output=True)
pfit2 = np.array(fit[0])

# Quartic-T, linear-C, quartic-N, linear-O:
guess = np.insert(pfit2, 7, [0, 0, 0])
fit = so.leastsq(error, guess, args=(logT,logC,logN,logO,logTOP,poly4141),
                 full_output=True)
pfit3 = np.array(fit[0])

# Quartic poly in all variables, fit only T<3000 K:
guess = list(pfit3[0: 6]) + [0,0,0] \
      + list(pfit3[6:11]) + [0,0,0]
fit = so.leastsq(error, guess, args=(logT,logC,logN,logO,logTOP,poly4, 3000),
                 full_output=True)
pfit4 = np.array(fit[0])
print("Polynomial coefficients for TOP:\n{}".format(pfit4))


# Show results:
Cmesh, Tmesh, Nmesh, Omesh = np.meshgrid(C, T, N, O)
model = 10**poly4(pfit4, logT, logC, logN, logO)

# Difference:
rand  = np.random.uniform(0.97,1.03,nt*nz**3)
rand2 = np.random.uniform(0.9,1.1,  nt*nz**3)
model[TOP==0] = 0.0
#model[Tmesh>5000] = TOP[Tmesh>5000]

# Percentage:
color = np.log(Nmesh.flatten())
plt.figure(0)
plt.clf()
ax = plt.subplot(111)
plt.xscale("log")
plt.scatter(Tmesh.flatten()*rand, 100*(1-model/TOP).flatten(), s=4.0, c=color)
plt.axhline( 10, ls="--", lw=1.5, zorder=-10, color="b")
plt.axhline(-10, ls="--", lw=1.5, zorder=-10, color="b")
plt.ylim(-100, 100)
plt.xlim(150, 7000)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([200, 300, 600, 1000, 2000, 3000, 6000])
plt.xlabel("Temperature  (K)")
plt.ylabel(r"(TOP-fit)/TOP  (%)")
plt.savefig("../plots/TOPpolyfit.pdf")
# Errors are typically less than 10% in pressure below 3000K
# Good enough for T < 3000 K
# points off at T~450K are blips in TEA

