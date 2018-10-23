#! /usr/bin/env python

import matplotlib.pyplot as plt


def text(x, y, text, fs):
  """
  Matplotlib text wrapper.
  """
  plt.text(x, y, text, fontsize=fs, ha="left", va="bottom",
           transform=ax.transAxes)

fs = 12

plt.figure(5, (6,8))
plt.clf()
plt.subplots_adjust(0.12, 0.25, 0.92, 0.92, wspace=0.17, hspace=0.5)
ax = plt.subplot(221)
plt.yscale("log")
plt.xlim(900, 3000)
plt.ylim(3e-2, 1e3)
rect = plt.Rectangle((2000, 10), 1000, 1000, fc="0.8", ec="k")
ax.add_patch(rect)
plt.xlabel("Temperature (K)", fontsize=fs)
plt.ylabel(r"$N_{\rm N}\,/\,N_{\rm C}$", fontsize=fs)
ax.set_xticks([2000])
ax.set_yticks([10])
ax.set_yticklabels(['10'])
ax.set_xticklabels(['2200'])
text(0.01, 1.01, r"$N_{\rm C}\,/\,N_{\rm O}< 0.1$", fs)
text(0.23, 0.32, r"HCO(CO)", fs)
text(0.55, 0.73, r"HCNO(CO)", fs)
ax.tick_params(labelsize=fs-1)

ax = plt.subplot(222)
plt.yscale("log")
plt.xlim(900, 3000)
plt.ylim(3e-2, 1e3)
rect = plt.Rectangle((2000, 10), 1000, 1000, fc="0.8", ec="k")
ax.add_patch(rect)
plt.xlabel("Temperature (K)", fontsize=fs)
ax.set_xticks([2000])
ax.set_yticks([10])
ax.set_yticklabels(['10'])
ax.set_xticklabels(['2200'])
text(0.01, 1.01, r"$0.1<N_{\rm C}\,/\,N_{\rm O}< 1$", fs)
text(0.23, 0.32, r"HCO(CO)", fs)
text(0.55, 0.73, r"HCNO(H$_2$O)", fs)
ax.tick_params(labelsize=fs-1)

ax = plt.subplot(223)
plt.yscale("log")
plt.xlim(200, 2300)
plt.ylim(1e-3, 3e2)
ax.add_patch(plt.Rectangle((900,  0.1), 2000, 1e3, fc="0.8", ec='k'))
plt.xlabel("Temperature (K)", fontsize=fs)
plt.ylabel(r"$N_{\rm N}\,/\,N_{\rm C}$", fontsize=fs)
ax.set_xticks([900])
ax.set_yticks([0.1])
ax.set_yticklabels(['0.1'])
text(0.01, 1.01, r"$N_{\rm C}\,/\,N_{\rm O}> 1,\qquad p < p_{\rm to}$", fs)
text(0.15, 0.18, r"HCO(H$_2$O)", fs)
text(0.45, 0.65, r"HCNO(H$_2$O)", fs)
ax.tick_params(labelsize=fs-1)

ax = plt.subplot(224)
plt.yscale("log")
plt.xlim(200, 2300)
plt.ylim(1e-3, 3e2)
plt.xlabel("Temperature (K)", fontsize=fs)
ax.set_xticks([900])
ax.set_yticks([0.1])
ax.set_yticklabels(['0.1'])
text(0.01, 1.01, r"$N_{\rm C}\,/\,N_{\rm O}> 1,\qquad p > p_{\rm to}$", fs)
text(0.35, 0.45, r"HCO(CO)", fs)
ax.tick_params(labelsize=fs-1)

plt.savefig("plots/regimes.ps")
