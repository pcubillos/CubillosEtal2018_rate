import sys
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../rate")
import rate

sys.path.append("../code")
import tea_wrapper as tw


def pprint(s):
  """
  Prettier print.
  """
  print(np.array2string(s, formatter={'int_kind':'{:6d}'.format}))


def plot(fign, C, N, O, c1, c2=None, c3=None, keep=False):
  """
  Plot abundances.
  """
  nlayers = len(c1[0])
  nmol    = len(c1)
  plt.figure(fign)
  if keep:
    for j in np.arange(nmol-1):
      plt.loglog(c1[j], press,   c=cols[j], lw=1.25, label=r.species[j])
      if c2 is not None:
        plt.loglog(c2[j], press, c=cols[j], lw=1.75, dashes=(10,2),zorder=-5)
      if c3 is not None:
        plt.loglog(c3[j], press, c="k", lw=2.0, dashes=(8,2,2,2),zorder=-10)
    return
  plt.clf()
  plt.loglog(np.tile(2*C, nlayers), press, dashes=(4,3), lw=2, color="0.5")
  plt.loglog(np.tile(2*N, nlayers), press, dashes=(4,3), lw=2, color="c")
  plt.loglog(np.tile(2*O, nlayers), press, dashes=(4,3), lw=2, color="r")
  for j in np.arange(nmol-1):
    plt.loglog(c1[j], press,   c=cols[j], lw=1.25, label=r.species[j])
    if c2 is not None:
      plt.loglog(c2[j], press, c=cols[j], lw=1.75, dashes=(10,2), zorder=-5)
    if c3 is not None:
      plt.loglog(c3[j], press, c="k", lw=2.0, dashes=(8,2,2,2), zorder=-10)
  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-20, 1)
  plt.legend(loc="upper left", framealpha=0.75, fontsize=10)
  plt.text(1e-19, 3e2, "T    = {:.0f} K\nC/O = {:.2g}\n"
           "N/C = {:.2g}\nM/H = {:.2g}".
           format(temp[0], C/O, N/C, C+N+O))
  plt.xlabel("Mixing ratio")
  plt.ylabel("Pressure (bar)")
  plt.axvspan(0.1, 1, facecolor='lightsteelblue', alpha=0.75, zorder=-201)


def qflag(sol, C, N, O, tol=1.01):
  """
  Check that abundances are within boundaries.
  """
  # Unpack sol:
  H2O, CH4, CO, CO2, NH3, C2H2, C2H4, HCN, N2, H2, H, He = sol
  f = (H + 2*H2) / H2

  flag = np.array([np.any(sol<0.0, axis=0),
                   (H2O /H2 > f*O  *tol),
                   (CH4 /H2 > f*C  *tol),
                   (CO  /H2 > f*C  *tol) | (CO /H2 > f*O  *tol),
                   (CO2 /H2 > f*C  *tol) | (CO2/H2 > f*O/2*tol),
                   (NH3 /H2 > f*N  *tol),
                   (C2H2/H2 > f*C/2*tol),
                   (C2H4/H2 > f*C/2*tol),
                   (HCN /H2 > f*C  *tol) | (HCN/H2 > f*N*tol),
                   (N2  /H2 > f*N/2*tol)])
  return flag


def main():
  # Setup:
  nlayers = 100
  press = np.logspace(-8, 3, nlayers)

  C, N, O = 2.5e-4, 1e-4, 5e-4
  r = rate.Rate(C=C, N=N, O=O, fHe=0.0)
  labels = r.species
  cols = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
          "pink", "0.5",    "gold",      "deepskyblue", "olive",   "seagreen"]

  nmol = len(r.species)
  nt   = 20
  nz   = 11
  temperature = np.logspace(2.305, 3.778, nt)
  carb = np.logspace(-6.57, -0.57, nz)
  nit  = np.logspace(-7.17, -1.17, nz)
  oxy  = np.logspace(-6.31, -0.31, nz)

  Cmesh, Tmesh, Nmesh, Omesh = np.meshgrid(carb, temperature, nit, oxy)


  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Explore solutions with C < O:
  flag1 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  flag2 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  flag3 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)

  hflag = np.zeros((nlayers, nt, nz,nz,nz), bool)

  print("This will take a few minutes.")
  ti = time.time()
  for t, T in enumerate(temperature):
    temp = np.tile(T, nlayers)
    for o, O in enumerate(oxy):
      for c, C in enumerate(carb):
        if C >= O:
          continue
        for n, N in enumerate(nit):
          if C + N + O > 0.1:
            continue
          sol1 = r.solve(temp, press, C, N, O, r.HCO_poly6_CO)
          sol2 = r.solve(temp, press, C, N, O, r.HCNO_poly8_CO)
          sol3 = r.solve(temp, press, C, N, O, r.HCNO_poly8_H2O)
          # Set flags:
          flag1[:,:,t,c,n,o] = qflag(sol1, C, N, O)
          flag2[:,:,t,c,n,o] = qflag(sol2, C, N, O)
          flag3[:,:,t,c,n,o] = qflag(sol3, C, N, O)
          hflag[:,t,c,n,o] = sol1[10] < sol1[9]  # H < H2
    print("{:2d}/{:d}".format(t+1, len(temperature)))
  tf = time.time()
  print(tf-ti)

  # Take a look at outputs, check when unphysical flags are raised, and
  # apply masks:
  print("Number of runs with unphysical abundances:\n"
  "Rows are HCO(CO), HCNO(CO), and HCNO(H2O), respectively.\n\n"
  "Without any masking:\n"
  "    Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2\n"
  +(":"*70))
  pprint(np.sum(flag1, axis=(1,2,3,4,5)))  # HCO CO
  pprint(np.sum(flag2, axis=(1,2,3,4,5)))  # HCNO CO
  pprint(np.sum(flag3, axis=(1,2,3,4,5)))  # HCNO H2O


  print("\nMask out high-temperature runs (where H2 dissociates into H):\n"
  "    Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2\n"
  +(":"*70))
  pprint(np.sum(flag1*hflag, axis=(1,2,3,4,5)))
  pprint(np.sum(flag2*hflag, axis=(1,2,3,4,5)))
  pprint(np.sum(flag3*hflag, axis=(1,2,3,4,5)))

  print("\nHCO: HCN flags as function of temperature (rows) and N (columns):")
  print(np.sum((flag1*hflag)[8], axis=(0,2,4)))
  print("At high temperatures, N/C>>1 and N/O>>1, HCO-HCN rise above\n"
        "upper C bounary at low pressures because poly does not account \n"
        "for HCN. Let's mask that out.\n")

  print("HCNO(CO): CO2 flags as function of temperature and N:")
  print(np.sum((flag2*hflag)[4], axis=(0,4,2)))
  print("HCNO(CO): CH4 flags as function of temperature and N:")
  print(np.sum((flag2*hflag)[2], axis=(0,4,2)))
  print("CO2 is overestimated at T~400-600 K.\n"
        "CH4 is overestimated at T<~700 K. Let's mask these out.\n")
  print("Abundances blip at T>~2000 K, high pressure, N/C>>1, and N/O>>1.\n"
        "Can't do much since HCO sol does not apply well in this regime.\n"
        "Not terribly worried since these are hardly expected conditions.\n")
  print("Actually, for some reason, HCNO(H2O) works better than HCNO(CO),\n"
        "when C/O>~0.1.  So, I'll take that solution in this regime.")

  Tlim = 2200.0
  NClim = 10.0
  COlim = 0.1
  print("\nApply HCO(CO), HCNO(CO), HCNO(H2O) at respecive ranges, delimited\n"
        "by boundaries in C/O ({:.1f}), N/C ({:.1f}), and T ({:.1f} K):\n".
          format(COlim,NClim,Tlim)
  +"    Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2\n"
  +(":"*70))
  pprint(np.sum(flag1*hflag*
                (1- (Nmesh/Cmesh > NClim)*(Tmesh > Tlim)), axis=(1,2,3,4,5)))
  pprint(np.sum(flag2*hflag * (Nmesh/Cmesh > NClim)
                            * (Cmesh/Omesh <= COlim)
                            * (Tmesh > Tlim), axis=(1,2,3,4,5)))
  pprint(np.sum(flag3*hflag * (Nmesh/Cmesh > NClim)
                            * (Cmesh/Omesh > COlim)
                            * (Tmesh > Tlim), axis=(1,2,3,4,5)))

  # Number of runs with unphysical abundances:
  # Rows are HCO(CO), HCNO(CO), and HCNO(H2O), respectively.

  # Without any masking:
  #     Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # [     0      0     10      0      0      0      0      0   2080      0]
  # [     0      0  96901      0  95158      0     10      0    141      0]
  # [ 52714      0 276554  46803   7731  18595  52613 142784     95  16890]
  #
  # Mask out high-temperature runs (where H2 dissociates into H):
  #     Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # [     0      0     10      0      0      0      0      0   1532      0]
  # [     0      0  96897      0   7842      0      8      0    139      0]
  # [ 52714      0 276554  46263   7731  18595  52613 142784     92  16890]
  #
  # Apply HCO(CO), HCNO(CO), HCNO(H2O) at respecive ranges, delimited
  # by boundaries in C/O (0.1), N/C (10.0), and T (2200.0 K):
  #     Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # [     0      0     10      0      0      0      0      0      8      0]
  # [     0      0    133      0      0      0      0      0     30      0]
  # [     0      0    277      0      0      0      0      0     15      0]


  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Explore solutions with C >= O:
  flag1 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  flag2 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  flag3 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  flag4 = np.zeros((10, nlayers, nt, nz,nz,nz), bool)
  hflag = np.zeros((nlayers, nt, nz,nz,nz), bool)

  print("\nBear with me again, a few more minutes.")
  ti = time.time()
  for t, T in enumerate(temperature):
    temp = np.tile(T, nlayers)
    for o, O in enumerate(oxy):
      for c, C in enumerate(carb):
        if C < O:
          continue
        for n, N in enumerate(nit):
          if C + N + O > 0.1:
            continue
          top = rate.top(temp, C, N, O)
          sol1 = r.solve(temp, press, C, N, O, r.HCO_poly6_CO)   # pnr
          sol2 = r.solve(temp, press, C, N, O, r.HCO_poly6_H2O)  # qnr
          sol3 = r.solve(temp, press, C, N, O, r.HCNO_poly8_CO)  # tnr
          sol4 = r.solve(temp, press, C, N, O, r.HCNO_poly8_H2O) # unr
          # Set flags:
          iCO = press > top
          flag1[:, iCO,t,c,n,o] = qflag(sol1[:, iCO],C,N,O)
          flag2[:,~iCO,t,c,n,o] = qflag(sol2[:,~iCO],C,N,O)
          flag3[:, iCO,t,c,n,o] = qflag(sol3[:, iCO],C,N,O)
          flag4[:,~iCO,t,c,n,o] = qflag(sol4[:,~iCO],C,N,O)
          hflag[:,t,c,n,o] = sol1[10] < sol1[9]  # H < H2
    print("{:2d}/{:d}".format(t+1, len(temperature)))
  tf = time.time()
  print(tf-ti)

  print("Number of runs with unphysical abundances:\n"
  "Rows are HCO(CO), HCNO(CO), HCO(H2O), and HCNO(H2O), respectively.\n\n"
  "Masking TOP, dissociating T:\n"
  "    Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2\n"
  +(":"*70))
  print("(Lower atmosphere, solve for CO:)")
  pprint(np.sum(flag1*hflag, axis=(1,2,3,4,5)))
  pprint(np.sum(flag3*hflag, axis=(1,2,3,4,5)))
  print("(Upper atmosphere, solve for H2O:)")
  pprint(np.sum(flag2*hflag, axis=(1,2,3,4,5)))
  pprint(np.sum(flag4*hflag, axis=(1,2,3,4,5)))

  # Masking TOP, dissociating T:
  #     Neg    H2O    CH4     CO    CO2    NH3   C2H2   C2H4    HCN     N2
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # (Lower atmosphere, solve for CO:)
  # [     0      0   1885      0      0      0      0      0    903      0]
  # [     0      0 179094      0      0      0  94992 108219  20232    227]
  # (Upper atmosphere, solve for H2O:)
  # [     0      0      0      0      0      0      0      0  19712      0]
  # [     0      0      0      0      0      0      0      0      0      0]


  print("CO:")
  print(np.sum(flag1[2], axis=(0,2,4)))
  print("HCO-CH4 fail at T~400-600 K and low pressure.\n"
        "But, values are barely above upper limit; so, don't worry.\n")
  print(np.sum(flag3[2], axis=(0,2,4)))
  print("HCNO-CH4 fails at low T, badly.\n")

  print("HCNO sol is too unstable;\n"
        "besides, CO/H2O>>1 disfavours HCN (bc low CH4),\n"
        "So, I'll simply take HCO solution.")


  print("\n\nH2O:\n"
  "I'm happy with both H2O solutions,\n"
  "HCO-HCN fails at T>~900 K and N/C>~0.1 (as expected),\n"
  "So, I'll take HCO when it works (bc it's faster), else take HCNO.")

  print("\nStill there are some blips here and there.\n"
        "Not much more than I can do.")


if __name__ == "__main__":
  main()

else:
  # Play around with some random isothermal profiles:
  # Copy/paste into an interactive Python session:

  # Setup:
  r = rate.Rate()
  nlayers = 100
  press = np.logspace(-8, 3, nlayers)
  labels = r.species
  cols = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
          "pink", "0.5",    "gold",      "deepskyblue", "olive",   "seagreen"]

  # Throw some random values for the abundances and temperature:
  O = 10**np.random.uniform(-6.31, -1.0)

  C = 10**np.random.uniform(-6.57, -1.0)
  #C = 10**np.random.uniform(-6.57, np.log10(O))  # C/O < 1.0
  #C = 10**np.random.uniform(np.log10(O), -0.57)  # C/O > 1.0

  N = 10**np.random.uniform(-7.17, -1.0)
  #N = 10**np.random.uniform(np.log10(10*C), -1.0)  # N/C > 10.0
  #N = 10**np.random.uniform(np.log10(0.1*C), -1.0)  # N/C > 0.1

  #temp = np.tile(10**np.random.uniform(2.305, 3.778), nlayers) # [200, 6000]
  temp = np.tile(10**np.random.uniform(2.305, 3.48), nlayers)  # [200, 3000]

  optimal  = r.solve(temp, press, C, N, O)
  HCO_CO   = r.solve(temp, press, C, N, O, r.HCO_poly6_CO)
  HCO_H2O  = r.solve(temp, press, C, N, O, r.HCO_poly6_H2O)
  HCNO_CO  = r.solve(temp, press, C, N, O, r.HCNO_poly8_CO)
  HCNO_H2O = r.solve(temp, press, C, N, O, r.HCNO_poly8_H2O)

  # This one should always look OK:
  plot(1, C, N, O, optimal)

  plot(2, C, N, O, HCO_CO)
  plt.title("HCO CO")

  plot(3, C, N, O, HCNO_CO)
  plt.title("HCNO CO")

  plot(4, C, N, O, HCO_H2O)
  plt.title("HCO H2O")

  plot(5, C, N, O, HCNO_H2O)
  plt.title("HCNO H2O")

  # Plot along with TEA:
  qtea = tw.tea(temp, press, C, N, O)
  plot(7, C,N,O, optimal, qtea)
