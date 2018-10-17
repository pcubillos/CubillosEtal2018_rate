#! /usr/bin/env python

import sys
import os
import time
import subprocess
import ctypes
import numpy as np
import scipy.constants as sc
import multiprocessing as mp

import tea_wrapper as tw


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
          qtea = tw.tea(T, press, C[c], N[n], O[o], t, c, n, o)
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
    temps = np.arange(i, nt, ncpu)
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
