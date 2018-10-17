#! /usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess


rootdir = os.path.realpath(os.path.dirname(__file__) + "/..")

def tea(temp, press, C, N, O, indices=None):
  """
  TEA wrapper to compute abundancecs for given atmospheric properties.

  Parameters
  ----------
  temp: 1D float iterable
     Atmospheric temperature profile (Kelvin).
  press: 1D float iterable
     Atmospheric pressure profile (bar).
  C: Float
     Carbon elemental fraction.
  N: Float
     Nitrogen elemental fraction.
  O: Float
     Oxygen elemental fraction.
  indices: 4-element integer tuple
     If not None, output-file indices (to differentiate runs).

  Returns
  -------
  q: list of 1D float arrays
     List of abundance profiles for H2O, CH4, CO, CO2, NH3, C2H2,
     C2H4, HCN, N2, H2, H.

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
  # Unpack indices for filenames:
  if indices is not None:
    t, c, n, o = indices
    patm = "preatm_{:02d}_{:02d}_{:02d}_{:02d}.tea".format(t,c,n,o)
    atmf = "atm_{:02d}_{:02d}_{:02d}_{:02d}".format(t,c,n,o)
  else:
    patm = "preatm.tea"
    atmf = "atm"

  # Make pre-atm:
  elements = ["H", "C", "N", "O"]
  nfrac    = [1.0,  C,   N,   O]
  specs = elements + ["H2", "H2O", "CH4", "CO", "CO2", "NH3",
                      "C2H2", "C2H4", "HCN", "N2"]
  d = np.loadtxt(rootdir+"/inputs/TEA_gdata_defaults.txt",
                 unpack=True, dtype=str)
  sd = {d[0,i]:d[1,i] for i in np.arange(len(d[0]))}
  sspecs = np.zeros(len(specs), "U50")
  for i, spec in enumerate(specs):
    sspecs[i] = sd[spec]

  with open(patm, 'w') as f:
    f.write("# TEA pre-atmosphere file.\n\n")
    f.write('#SPECIES\n{:s}\n\n'.format(" ".join(sspecs)))
    f.write("#TEADATA\n#Pressure          Temp  " +
          "  ".join(["{:>12s}".format(atom) for atom in elements])+"\n")
    for p,t in zip(press,temp):
      f.write("{:10.4e}     {:>8.2f}  ".format(p,t))
      f.write("  ".join(["{:12.6e}".format(abun) for abun in nfrac]) + "\n")
  # Run TEA:
  proc = subprocess.Popen([rootdir+"/TEA/tea/runatm.py", patm, atmf])
  proc.communicate()
  # Gather results:
  d = np.loadtxt(atmf+".tea", skiprows=8, unpack=True)
  os.remove(patm)
  #      H2O,  CH4,  CO,   CO2,   NH3,   C2H2,  C2H4,  HCN,   N2,    H2,   H:
  return d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15], d[6], d[2]
