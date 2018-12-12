#! /usr/bin/env python

import sys
import ctypes
import numpy as np
import multiprocessing as mp

import tea_wrapper as tw

sys.path.append("../rate")
import rate


def worker(H2O, CO, temp, press, C, N, O, itemps):
  """
  Multiprocessing wrapper to run TEA calculations.
  """
  # Rate object to provide initial guesses for TEA:
  r = rate.Rate(C=1e-4, N=1e-4, O=1e-4, fHe=0.0)

  for c in range(len(C)):
    for o in range(len(O)):
      if C[c] < O[o]:
        continue
      for n in range(len(N)):
        for t in itemps:
          T = np.tile(temp[t], len(press))
          tail = "_{:02d}_{:02d}_{:02d}_{:02d}".format(t,c,n,o)
          sol = r.solve(T, press, C[c], N[n], O[o])
          guess = {mol: val for mol, val in zip(r.species, sol[:,0])}
          qtea = tw.tea(T, press, C[c], N[n], O[o], tail, guess)
          H2O[:,t,c,n,o] = qtea[0]
          CO [:,t,c,n,o] = qtea[2]
          if 0 in itemps:
            if len(itemps) == 1:
              dt = 1
            else:
              dt = itemps[1] - itemps[0]
            print("[{:2d}/{:d} {:2d}/11 {:2d}/11 {:2d}/11]"
                  .format(t//dt, len(itemps), c, o, n))


def main():
  # Number of layers, temperatures, and elemental abundances:
  nl = 100
  nt = 20
  nz = 11
  # Number of parallel CPUs:
  ncpu = 8

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
