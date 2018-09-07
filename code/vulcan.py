#! /usr/bin/env python

import sys
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../vulcan")
import vulcan_analytic as v

col = ["orange", "dodgerblue", "m", "c", "g", "k", "y", "gold", "r", "b", "navy"]
allspec = ["H2", "He", "H2O", "CO", "CO2", "CH4", "C2H2", "C2H4", "HCN", "NH3", "N2"]

def qplot(q, fignum, yran=[1e-25, 1], savefile=None):
  plt.figure(fignum)
  plt.clf()
  for i in np.arange(len(allspec)):
    plt.semilogy(temp, q[:,i], lw=2, color=col[i], label=allspec[i])
  plt.ylim(yran)
  plt.xlim(500, 3000)
  plt.legend(loc="lower right", fontsize=9)
  plt.xlabel("Temperature (K)")
  plt.ylabel("Mixing fraction")
  if savefile is not None:
    plt.savefig(savefile)


nlayers = 100
temp = np.linspace(500.0, 3000.0, nlayers)
press = np.tile(1.0, nlayers)

# Fig 1, top panel of Heng & Tsai (2016) works fine:
ti = time.time()
q = v.run_analytic_vulcan(temp, press, n_o=5e-4, n_c=2.5e-4, n_n=1e-4)
print("VULCAN 100 layers: {:.4f} s.".format(time.time()-ti))
qplot(q, 1)

# Fig 1, middle panel of Heng & Tsai (2016) works fine too:
ti = time.time()
q = v.run_analytic_vulcan(temp, press, n_o=5e-4, n_c=5e-4, n_n=1e-4)
print("VULCAN 100 layers: {:.4f} s.".format(time.time()-ti))
qplot(q, 2)

# Lowering the pressure starts to cause artifacts:
press = np.tile(0.01, nlayers)
ti = time.time()
q = v.run_analytic_vulcan(temp, press, n_o=5e-4, n_c=5e-4, n_n=1e-4)
print("VULCAN 100 layers: {:.4f} s.".format(time.time()-ti))
qplot(q, 3, yran=[1e-30,1],
      savefile="../plots/VULCAN_6mol-poly_0.01bar_O5e-4_C5e-4_N1e-4.pdf")

# Artifacts increase with lower pressure:
press = np.tile(1e-4, nlayers)
ti = time.time()
q = v.run_analytic_vulcan(temp, press, n_o=5e-4, n_c=5e-4, n_n=1e-4)
print("VULCAN 100 layers: {:.4f} s.".format(time.time()-ti))
qplot(q, 4, yran=[1e-30,1])

