#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../code")
import tea_wrapper as tw

sys.path.append("../rate")
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
qtea1 = tw.tea(temp, press, carbon1, nitrogen, oxygen)
# C/O > 1.0:
Q2 = r.solve(temp, press, C=carbon2)
qtea2 = tw.tea(temp, press, carbon2, nitrogen, oxygen)

fs = 14
thin = 4
xran = 1e-20, 2.0
labels = r.species
cols   = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
          "pink", "0.5",    "gold",      "deepskyblue", "olive",    "seagreen"]


plt.figure(-1, (6,8))
plt.clf()
# C/O < 1:
ax = plt.axes([0.12, 0.57, 0.55, 0.4])
for q, qtea, col, lab in zip(Q1, qtea1, cols, labels):
  plt.loglog(q, press, color=col, lw=2.0, label=lab, zorder=-2)
  plt.loglog(qtea[::thin], press[::thin], ls="", lw=0.5, marker="o",
             ms=3.5, color=col, mec="k", mew=0.5)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.legend(loc="upper left", framealpha=0.75, fontsize=fs-5)
plt.text(1.5e-20, 4e2, "C/O = {:.1f}\n".format(carbon1/oxygen), fontsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
ax.tick_params(labelsize=fs-2)
ax.set_xticks([1e-20, 1e-15, 1e-10, 1e-5, 1])
# Residuals:
ax = plt.axes([0.69, 0.57, 0.29, 0.4])
for q, qtea, col in zip(Q1, qtea1, cols):
  plt.loglog(100*np.abs(1-q/qtea), press, lw=1.5, color=col)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-5, 1e2)
ax.set_yticklabels([])
ax.tick_params(labelsize=fs-2)
ax.set_xticks([1e-4, 1e-2, 1, 100])
plt.xlabel("Difference (%)", fontsize=fs)

# C/O > 1:
ax = plt.axes([0.12, 0.1, 0.55, 0.4])
for q, qtea, col, lab in zip(Q2, qtea2, cols, labels):
  plt.loglog(q, press, color=col, lw=2.0, label=lab, zorder=-2)
  plt.loglog(qtea[::thin], press[::thin], ls="", lw=0.5, marker="o",
             ms=3.5, color=col, mec="k", mew=0.5)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.text(1.5e-20, 4e2, "C/O = {:.1f}\n".format(carbon2/oxygen), fontsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)
plt.ylabel("Pressure (bar)", fontsize=fs)
ax.tick_params(labelsize=fs-2)
ax.set_xticks([1e-20, 1e-15, 1e-10, 1e-5, 1])
# Residuals:
ax = plt.axes([0.69, 0.1, 0.29, 0.4])
for q, qtea, col in zip(Q2, qtea2, cols):
  plt.loglog(100*np.abs(1-q/qtea), press, lw=1.5, color=col)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-5, 1e2)
ax.set_yticklabels([])
ax.set_xticks([1e-4, 1e-2, 1, 100])
ax.tick_params(labelsize=fs-2)
plt.xlabel("Difference (%)", fontsize=fs)

plt.savefig("../plots/bench.ps")

