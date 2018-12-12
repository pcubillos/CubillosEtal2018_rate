import sys
import numpy as np
import multiprocessing as mp

sys.path.append("../rate")
import rate

sys.path.append("../code")
import tea_wrapper as tw


def worker(temps, press, C, N, O, MH, C2O, itemps):
  """
  Multiprocessing wrapper to run TEA calculations.
  """
  nmetal, nco = np.shape(C)
  # Rate object to provide initial guesses for TEA:
  r = rate.Rate(C=1e-4, N=1e-4, O=1e-4, fHe=0.0)

  for iz in range(nmetal):
    for ico in range(nco):
      for it in itemps:
        T = np.tile(temps[it], len(press))
        tail = "_metal{:+.1f}_COrat{:.2f}_temp{:02d}".format(
                   MH[iz], C2O[ico], it)
        sol = r.solve(T, press, C[iz,ico], N[iz,ico], O[iz,ico])
        guess = {mol: val for mol, val in zip(r.species, sol[:,0])}
        qtea = tw.tea(T, press, C[iz,ico], N[iz,ico], O[iz,ico], tail, guess)
        if 0 in itemps:
          print("[{:2d}/{:d}  {:2d}/{:d}  {:2d}/{:d}]"
                .format(it, len(itemps), iz, nmetal, ico, nco))


def main():
  # Solar values:
  Csun = 2.69e-4
  Nsun = 6.76e-5
  Osun = 4.90e-4
  Zsun = Csun + Nsun + Osun
  # Make grid of constant C/O vs constant C+N+O:
  MH  = np.array([-3.0, 0.0, 1.0, 2.0, 3.0])
  Z   = 10**MH * Zsun
  C2O = np.array([0.1, 0.55, 0.9, 5.0])
  nmetal = len(MH)
  nco    = len(C2O)

  # Elemental abundances for this grid:
  C = np.zeros((nmetal,nco))
  N = np.zeros((nmetal,nco))
  O = np.zeros((nmetal,nco))
  for iz,z in enumerate(Z):
    for ico,c2o in enumerate(C2O):
      C[iz,ico] = z / (1.0 + 1.0/c2o + 1.0/(Csun/Nsun))
      N[iz,ico] = C[iz,ico] / (Csun/Nsun)
      O[iz,ico] = C[iz,ico] / c2o

  # TEA grid:
  nl = 100
  nt = 20
  # Number of parallel CPUs:
  ncpu = 8

  press = np.logspace(-8, 3, nl)
  temps = np.logspace(2.305, 3.778, nt)

  procs = []
  for i in range(ncpu):
    itemps = np.arange(i, nt, ncpu)
    p = mp.Process(target=worker,
                   args=(temps, press, C, N, O, MH, C2O, itemps))
    p.start()
    procs.append(p)

  # Wait until they compute the whole grid:
  for i in np.arange(ncpu):
    procs[i].join()


if __name__ == "__main__":
  main()
