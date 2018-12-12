import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../rate")
import rate


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make grid of constant C/O vs constant C+N+O (be consistent with benchmark.py):
Csun = 2.69e-4
Nsun = 6.76e-5
Osun = 4.90e-4
Zsun = Csun + Nsun + Osun

# Define grid:
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

# Atmosphere grid:
nlayers = 100
nt = 20
press = np.logspace(-8, 3, nlayers)

He = np.zeros(nlayers)

press = np.logspace(-8, 3, nlayers)
temps = np.logspace(2.305, 3.778, nt)

r = rate.Rate(C=2.5e-4, N=1e-4, O=5e-4, fHe=0.0)


# Gather TEA results and evaluate rate over grid:
tea_grid  = np.zeros((len(r.species), nlayers, nt, nmetal, nco))
rate_grid = np.zeros((len(r.species), nlayers, nt, nmetal, nco))

for it, T in enumerate(temps):
  temp = np.tile(T, nlayers)
  for iz in range(nmetal):
    for ico in range(nco):
      tail = "_metal{:+.1f}_COrat{:.2f}_temp{:02d}".format(MH[iz], C2O[ico], it)
      d = np.loadtxt("atm{:s}.tea".format(tail), skiprows=8, unpack=True)
      #                    H2O, CH4, CO, CO2, NH3, C2H2, C2H4, HCN, N2, H2, H:
      tea_grid [:,:,it,iz,ico] = [d[ 7], d[ 8], d[ 9], d[10], d[11], d[12],
                                  d[13], d[14], d[15], d[ 6], d[ 2], He]
      rate_grid[:,:,it,iz,ico] = r.solve(temp, press, C[iz,ico], N[iz,ico],
                                         O[iz,ico])


# The plot:
TT, PP = np.meshgrid(temps, press)
xticks = [300, 1000, 5000]

logdiff = np.abs(np.log10(rate_grid/tea_grid))

fs = 8.25
labels = [r'H$_2$O',   r'CH$_4$', r'CO',    r'CO$_2$', r'NH$_3$', r'C$_2$H$_2$',
        r'C$_2$H$_4$', r'HCN',    r'N$_2$', r'H$_2$',  r'H',      r'He']
cols = ["navy", "orange", "limegreen", "red",         "magenta", "brown",
        "pink", "0.5",    "gold",      "deepskyblue", "olive",   "seagreen"]

# logdiff contour levels:
#       0,   10%     50%    100%  10x  1000x
levs = [0.0, 0.043,  0.176, 0.3,  1.0, 3.0]
cols = plt.cm.viridis_r(np.linspace(0, 0.9, len(levs)))
cols = plt.cm.viridis_r(np.array([0.  , 0.18, 0.36, 0.54, 0.72, 0.9,0.9 ]))

for imol in range(len(r.species)-1):
  plt.figure(-4, (8.5, 5.4))

  plt.clf()
  plt.subplots_adjust(0.08, 0.09, 0.88, 0.95, wspace=0.1, hspace=0.1)
  for iz,z in enumerate(MH):
    for ic,c in enumerate(C2O):
      ax = plt.subplot(nmetal, nco, ic+1+nco*iz)
      contour = np.clip(logdiff[imol,:,:,iz,ic], 0, 3)
      csf = ax.contourf(TT, PP, contour, zorder=1, levels=levs[:],
                        colors=cols, extend='max')
      hatch = tea_grid[imol,:,:,iz,ic]
      hcs = plt.contourf(TT, PP, hatch, levels=[0, 1e-10], colors='none',
                        hatches=['\\\\', None], extend='lower')
      ax.set_xscale('log')
      ax.set_yscale('log')
      ax.set_ylim(np.amax(press), np.amin(press))
      ax.set_yticks([1e2, 1e-3, 1e-8])
      ax.set_xticks(xticks)
      ax.set_xticklabels(xticks)
      ax.tick_params(labelsize=fs)
      if ic > 0:
        ax.set_yticklabels([])
      else:
        plt.ylabel("Pressure (bar)", fontsize=fs)
      if iz == 0:
        ax.set_title("C/O = ${:.2g}$".format(c), fontsize=fs)
      if ic == nco-1:
        ax.text(1.05, 0.5, "[M/H] = ${:.1f}$".format(z), fontsize=fs,
                transform=ax.transAxes, va='center', rotation='vertical')
      if iz == 1:
        plt.text(1.195, 0.5, "{:s}".format(labels[imol]), fontsize=2.5*fs,
                 transform=ax.transAxes, va="bottom")
      if ic == nco-1 and iz==0:
        bounds = levs + [3.1]
        pos = ax.get_position()
        (x0, y0),(x1,y1) = ax.get_position().get_points()
        cb = plt.colorbar(csf, cax=plt.axes([x1+0.04, y0, 0.01, y1-y0]),
            extend='max', boundaries=levs, ticks=levs[:-1])
        cb.set_ticklabels(['0', '10%', '50%', '100%', r'10$\times$', ''])
      if iz == nmetal-1:
        plt.xlabel("Temperature (K)", fontsize=fs)
      else:
        ax.set_xticklabels([])

  plt.savefig("../plots/benchmark_rate-tea_{:s}.pdf".format(r.species[imol]))
  plt.savefig("../plots/benchmark_rate-tea_{:s}.ps".format(r.species[imol]))
