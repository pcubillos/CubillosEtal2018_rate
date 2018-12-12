#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("rate")
import rate as rate


nlayers = 100
temp  = np.tile(1200.0, nlayers)
press = np.logspace(-8, 3, nlayers)

carbon1  = 2.5e-4
carbon2  = 1.0e-3
nitrogen = 7.0e-5
oxygen   = 5.0e-4
r = rate.Rate(C=carbon1, N=nitrogen, O=oxygen)
# C/O < 1.0:
Q1 = r.solve(temp, press)
# C/O > 1.0:
Q2 = r.solve(temp, press, C=carbon2)
top = rate.top(temp[0], carbon2, nitrogen, oxygen)

fs = 14
xran = 1e-12, 1e-2
labels = r.species
cols  = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
         "pink", "0.5",    "gold",      "deepskyblue", "olive",   "seagreen"]

plt.figure(-1, (6,8))
plt.clf()
plt.subplots_adjust(0.13, 0.08, 0.95, 0.97, hspace=0.25)
ax = plt.subplot(211)
plt.loglog(np.tile(2*carbon1,nlayers), press, dashes=(4,2), lw=1.75,
           color="k", zorder=-3)
plt.loglog(np.tile(2*oxygen,nlayers), press, dashes=(4,2), lw=1.75,
           color="r", zorder=-3)
for q, lab, col in zip(Q1, labels, cols):
  if lab in ["H", "H2", "He"]:
    continue
  elif lab in ["H2O", "CO"]:
    lw = 2.0
    dashes = ()
  else:
    lw = 1.5
    dashes = (10,1,3,1)
  plt.loglog(q, press, color=col, lw=lw, dashes=dashes, label=lab)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.legend(loc="upper left", framealpha=0.75, fontsize=fs-4)
plt.text(1.5e-10, 5e-7, "C/O = {:.1f}\n".format(carbon1/oxygen), fontsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
ax.tick_params(labelsize=fs-2)

ax = plt.subplot(212)
plt.loglog(np.tile(2*carbon2,nlayers), press, dashes=(4,2), lw=1.75,
           color="k", zorder=-3)
plt.loglog(np.tile(2*oxygen,nlayers), press, dashes=(4,2), lw=1.75,
           color="r", zorder=-3)
for q, lab, col in zip(Q2, labels, cols):
  if lab in ["H", "H2", "He"]:
    continue
  elif lab in ["H2O", "CO"]:
    lw = 2.0
    dashes = ()
  else:
    lw = 1.5
    dashes = (10,1,3,1)
  plt.loglog(q, press, color=col, lw=lw, dashes=dashes, label=lab)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.text(2e-12, 5e-7, "C/O = {:.1f}\n".format(carbon2/oxygen), fontsize=fs-2)
plt.xlabel("Mixing ratio",   fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
plt.axhline(top, ls="-", lw=4, color="0.7", zorder=-101)
plt.text(1.2e-12,  1, "CO dominated",  fontsize=fs-4)
plt.text(1.2e-12, 10, "H2O dominated", fontsize=fs-4)
ax.tick_params(labelsize=fs-2)

plt.savefig("plots/domains.ps")
plt.savefig("plots/domains.pdf")
