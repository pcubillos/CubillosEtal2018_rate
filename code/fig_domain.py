#! /usr/bin/env python

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sys.path.append("rate")
import rate as rate


#C, N, O = 2.5e-4, 7.0e-5, 5.0e-4
#r = rate.Rate(C, N, O)

r = rate.Rate()
C = r.C = 2.5e-4
N = r.N = 7e-5
O = r.O = 5e-4

nlayers = 100
temp  = np.tile(1200.0, nlayers)
press = np.logspace(-8, 3, nlayers)

fs = 14
xran = 1e-12, 1e-2
label = "H2O",  "CH4",  "CO",     "CO2", "NH3", "C2H2", "C2H4", "HCN", "N2"
col   = "navy", "orange", "limegreen", "r", "m", "brown", "pink", "0.5", "gold"

plt.figure(-1, (6,8))
plt.clf()
plt.subplots_adjust(0.13, 0.08, 0.95, 0.97, hspace=0.25)
# C/O < 1.0:
Q1 = r.solve(temp, press, C=2.5e-4)
ax = plt.subplot(211)
plt.loglog(np.tile(2*r.C,nlayers), press, dashes=(4,2), lw=1.75,
           color="k", zorder=-3)
plt.loglog(np.tile(2*r.O,nlayers), press, dashes=(4,2), lw=1.75,
           color="r", zorder=-3)
plt.loglog(Q1[0], press,   color=col[0], lw=2.0, label=label[0])
plt.loglog(Q1[2], press,   color=col[2], lw=2.0, label=label[2])
for j in np.arange(len(Q1)):
  if j not in [0,2]:
    plt.loglog(Q1[j], press, color=col[j], lw=1.25, label=label[j],
               dashes=(10,1,3,1), zorder=-2)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.legend(loc="upper left", framealpha=0.75, fontsize=10)
plt.text(1.5e-10, 5e-7, "C/O = {:.2g}\n".format(r.C/r.O), fontsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
ax.tick_params(labelsize=fs-2)

# C/O > 1.0:
Q2 = r.solve(temp, press, C=1e-3)
top = rate.top(temp[0], r.C, r.N, r.O)
ax = plt.subplot(212)
plt.loglog(np.tile(2*r.C,nlayers), press, dashes=(4,2), lw=1.75,
           color="k", zorder=-3)
plt.loglog(np.tile(2*r.O,nlayers), press, dashes=(4,2), lw=1.75,
           color="r", zorder=-3)
plt.loglog(Q2[0], press,   color=col[0], lw=2.0, label=label[0])
plt.loglog(Q2[2], press,   color=col[2], lw=2.0, label=label[2])
for j in np.arange(len(Q2)):
  if j not in [0,2]:
    plt.loglog(Q2[j], press, color=col[j], lw=1.25, label=label[j],
               dashes=(10,1,3,1), zorder=-2)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.text(2e-12, 5e-7, "C/O = {:.2g}\n".format(r.C/r.O), fontsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
plt.axhline(top, ls="-", lw=4, color="0.7", zorder=-101)
plt.text(1.2e-12,  1, "CO dominated", fontsize=fs-4)
plt.text(1.2e-12, 10, "H2O dominated", fontsize=fs-4)
ax.tick_params(labelsize=fs-2)

plt.savefig("plots/domains.ps")
