#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../code")
import tea_wrapper as tw

sys.path.append("../rate")
import rate as rate


nlayers = 100
temp1 = np.tile(1900.0, nlayers)
temp2 = np.tile(2500.0, nlayers)
temp3 = np.tile(3000.0, nlayers)
press = np.logspace(-8, 3, nlayers)

carbon   = 1.5e-4
nitrogen = 7.0e-5
oxygen   = 5.0e-4
r = rate.Rate(C=carbon, N=nitrogen, O=oxygen, fHe=0.0)
Q1 = r.solve(temp1, press)
Q2 = r.solve(temp2, press)
Q3 = r.solve(temp3, press)
qtea1 = tw.tea(temp1, press, carbon, nitrogen, oxygen)
qtea2 = tw.tea(temp2, press, carbon, nitrogen, oxygen)
qtea3 = tw.tea(temp3, press, carbon, nitrogen, oxygen)


fs = 11
xran = 1e-15, 2

labels = r.species
cols = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
        "pink", "0.5",    "gold",      "deepskyblue", "olive",   "seagreen"]

ishow = np.where(np.in1d(labels, ['H2O', 'CO', 'CO2', 'HCN', 'H2', 'H']))[0]
xtext = 0.03

plt.figure(-20, (9, 4.5))
plt.clf()
plt.subplots_adjust(0.1, 0.45, 0.95, 0.97, hspace=0.25, wspace=0.23)
ax = plt.subplot(131)
for j in ishow:
  plt.loglog(Q1[j],    press, c=cols[j], lw=1.5, label=labels[j])
  plt.loglog(qtea1[j], press, c=cols[j], lw=1.0, dashes=(6,2,3,2))

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.legend(framealpha=0.75, fontsize=fs-4, bbox_to_anchor=(0.3, 0.55))
plt.text(xtext, 0.03, r"$T=1900\,$K", fontsize=fs-2, transform=ax.transAxes)
plt.ylabel("Pressure (bar)", fontsize=fs)
ax.tick_params(labelsize=fs-2)
plt.xlabel("Mixing ratio", fontsize=fs)

ax = plt.subplot(132)
for j in ishow:
  plt.loglog(Q2[j], press, color=cols[j], lw=1.5, label=labels[j])
  plt.loglog(qtea2[j], press, c=cols[j], lw=1.0, dashes=(6,2,3,2))
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.text(xtext, 0.03, r"$T=2500\,$K", fontsize=fs-2, transform=ax.transAxes)
plt.xlabel("Mixing ratio", fontsize=fs)

ax = plt.subplot(133)
for j in ishow:
  plt.loglog(Q3[j], press, color=cols[j], lw=1.5, label=labels[j])
  plt.loglog(qtea3[j], press, c=cols[j], lw=1.0, dashes=(6,2,3,2))
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(xran)
plt.text(xtext, 0.03, r"$T=3000\,$K", fontsize=fs-2, transform=ax.transAxes)
plt.xlabel("Mixing ratio", fontsize=fs)

plt.savefig("../plots/hot_limit.ps")
