#! /usr/bin/env python

import os
import subprocess
import numpy as np


rootdir = os.path.realpath(os.path.dirname(__file__) + "/..")

def tea(temp, press, C, N, O, tail='', guess={}):
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
  tail: string
     output file name tail.
  guess: dict
     Dictionary with guess values for species.

  Returns
  -------
  q: list of 1D float arrays
     List of abundance profiles for H2O, CH4, CO, CO2, NH3, C2H2,
     C2H4, HCN, N2, H2, H.

  Example
  -------
  >>> import tea_wrapper as tw
  >>> nlayers = 100
  >>> press = np.logspace(-8, 3, nlayers)
  >>> temp  = np.tile(1500.0, nlayers)
  >>> C, N, O = 2.5e-4, 1e-4, 5e-4
  >>> qtea = tw.tea(temp, press, C, N, O)
  """
  # File names:
  patm   = "preatm{:s}.tea".format(tail)
  atmf   = "atm{:s}".format(tail)
  guessf = "guess{:s}.txt".format(tail)

  # Make pre-atm:
  elements = ["H", "C", "N", "O"]
  nfrac    = [1.0,  C,   N,   O]
  specs = elements + ["H2", "H2O", "CH4", "CO", "CO2", "NH3",
                      "C2H2", "C2H4", "HCN", "N2"]
  d = np.array(np.loadtxt(rootdir+"/inputs/TEA_gdata_defaults.txt",
                          unpack=True, dtype=bytes), "U50")
  sd = {d[0,i]:d[1,i] for i in np.arange(len(d[0]))}
  sspecs = np.zeros(len(specs), "U50")
  guesses = np.tile(1e-50, len(specs))
  for i, spec in enumerate(specs):
    sspecs[i] = sd[spec]
    if spec in guess:
      guesses[i] = guess[spec]

  with open(guessf, 'w') as f:
    f.write(" ".join("{:.5e}".format(val) for val in guesses))

  with open(patm, 'w') as f:
    f.write("# TEA pre-atmosphere file.\n\n")
    f.write('#SPECIES\n{:s}\n\n'.format(" ".join(sspecs)))
    f.write("#TEADATA\n#Pressure          Temp  " +
          "  ".join(["{:>12s}".format(atom) for atom in elements])+"\n")
    for p,t in zip(press,temp):
      f.write("{:10.4e}     {:>8.2f}  ".format(p,t))
      f.write("  ".join(["{:12.6e}".format(abun) for abun in nfrac]) + "\n")

  # Run TEA:
  if guess != {}:
    proc = subprocess.Popen([rootdir+"/TEA/tea/runatm.py", patm, atmf, guessf])
  else:
    proc = subprocess.Popen([rootdir+"/TEA/tea/runatm.py", patm, atmf])
  proc.communicate()
  # Gather results:
  d = np.loadtxt(atmf+".tea", skiprows=8, unpack=True)
  os.remove(patm)
  os.remove(guessf)
  #      H2O,  CH4,  CO,   CO2,   NH3,   C2H2,  C2H4,  HCN,   N2,    H2,   H:
  return d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15], d[6], d[2]
