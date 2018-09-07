#! /usr/bin/env python

import sys
import os
import time
import subprocess
import ctypes
import numpy as np
import scipy.constants as sc
import multiprocessing as mp


def tea(temp, press, C, N, O, t, c, n, o):
  """
  TEA wrapper to compute abundancecs for given atmospheric properties.

  Example
  -------
  >>> nlayers = 100
  >>> press = np.logspace(-8, 3, nlayers)
  >>> temp  = np.tile(1500.0, nlayers)
  >>> C, N, O = 2.5e-4, 1e-4, 5e-4

  >>> ti = time.time()
  >>> qtea = tea(temp, press, C, N, O)
  >>> print(time.time()-ti)
  """
  # Make pre-atm:
  elements = ["H", "C", "N", "O"]
  nfrac    = [1.0,  C,   N,   O]
  specs = elements + ["H2", "H2O", "CH4", "CO", "CO2", "NH3",
                      "C2H2", "C2H4", "HCN", "N2"]
  d = np.loadtxt("../inputs/TEA_gdata_defaults.txt", unpack=True, dtype=str)
  sd = {d[0,i]:d[1,i] for i in np.arange(len(d[0]))}
  sspecs = np.zeros(len(specs), "|S50")
  for i in np.arange(len(specs)):
    sspecs[i] = sd[specs[i]]
  patm = "./preatm_{:02d}_{:02d}_{:02d}_{:02d}.tea".format(t,c,n,o)
  with open(patm, 'w') as f:
    f.write("# TEA pre-atmosphere file.\n\n")
    f.write('#SPECIES\n{:s}\n\n'.format(" ".join(sspecs)))
    f.write("#TEADATA\n#Pressure          Temp  " +
          "  ".join(["{:>12s}".format(atom) for atom in elements])+"\n")
    for i in np.arange(len(temp)):
      f.write("{:10.4e}     {:>8.2f}  ".format(press[i], temp[i]))
      f.write("  ".join(["{:12.6e}".format(abun) for abun in nfrac]) + "\n")
  # Run TEA:
  atmf = "atm_{:02d}_{:02d}_{:02d}_{:02d}".format(t,c,n,o)
  proc = subprocess.Popen(["../TEA/tea/runatm.py", patm, atmf])
  proc.communicate()
  # Gather results:
  d = np.loadtxt(atmf+".tea", skiprows=8, unpack=True)
  os.remove(patm)
  #      H2O,  CH4,  CO,   CO2,   NH3,   C2H2,  C2H4,  HCN,   N2,    H2:
  return d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15], d[6]


def worker(H2O, CO, temp, press, C, N, O, temps):
  """
  Multiprocessing wrapper to run TEA calculations.
  """
  for c in np.arange(len(C)):
    for o in np.arange(len(O)):
      if C[c] < O[o]:
        continue
      for n in np.arange(len(N)):
        for t in temps:
          T = np.tile(temp[t], len(press))
          qtea = tea(T, press, C[c], N[n], O[o], t, c, n, o)
          H2O[:,t,c,n,o] = qtea[0]
          CO [:,t,c,n,o] = qtea[2]
          if 0 in temps:
            print("[{:2d}/20 {:2d}/11 {:2d}/11 {:2d}/11]".format(t, c, o, n))


def main():
  # Number of layers, temperatures, and elemental abundances:
  nl = 100
  nt = 20
  nz = 11
  # Number of parallel CPUs:
  ncpu = 7

  press = np.logspace(-8, 3, nl)
  temp = np.logspace(2.305, 3.778, nt)
  C    = np.logspace(-6.57, -0.57, nz)
  N    = np.logspace(-7.17, -1.17, nz)
  O    = np.logspace(-6.31, -0.31, nz)

  # Shared memory arrays:
  sm_H2O = mp.Array(ctypes.c_double, nl*nt*nz*nz*nz)
  sm_CO  = mp.Array(ctypes.c_double, nl*nt*nz*nz*nz)
  H2O = np.ctypeslib.as_array(sm_H2O.get_obj()).reshape((nl,nt,nz,nz,nz))
  CO  = np.ctypeslib.as_array(sm_CO.get_obj()).reshape( (nl,nt,nz,nz,nz))

  procs = []
  for i in np.arange(ncpu):
    temps = np.arange(i,nt, ncpu)
    p = mp.Process(target=worker, args=(H2O, CO, temp, press, C, N, O, temps))
    p.start()
    procs.append(p)

  # Wait until they compute the whole grid:
  for i in np.arange(ncpu):
    procs[i].join()

  # Save to npz file:
  np.savez("TEA_grid.npz", H2O=H2O, CO=CO)


if __name__ == "__main__":
  main()
