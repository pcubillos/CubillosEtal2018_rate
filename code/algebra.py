#! /usr/bin/env python

import numpy as np
import sympy as s


def get_coeffs(eq, var):
  """
  Get algebraic formula for the polynomial coefficients of an equation.
  """
  eq = s.simplify(eq)
  numer = s.numer(eq)
  denom = s.denom(eq)
  eqq = numer*denom

  p = s.Poly(numer, var)
  return p.coeffs()

# Define algebraic variable:
k1, k2, k3, k4, k5, k6, f = s.symbols("k1 k2 k3 k4 k5 k6 f")
C, N, O = s.symbols("C N O")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Heng & Lyon (2016):
CH4 = s.symbols("CH4")

C2H2 = k3 * CH4**2  # Eq 0

# Mass balance:
# CH4 + CO + CO2 + 2*C2H2 = 2*C  # Eq1
# H2O + CO + 2*CO2        = 2*O  # Eq2

## Solve for H2O from:  2*Eq1 - Eq2
#    2*CH4 + CO         + 4*C2H2 - H2O = 4*C - 2*O
# => 2*CH4 + k1*CH4*H2O + 4*C2H2 - H2O = 4*C - 2*O
# => H2O * (k1*CH4-1) = 4*C - 2*O - 2*CH4 - 4*C2H2
H2O = (4*C - 2*O - 2*CH4 - 4*C2H2) / (k1*CH4-1)
CO  = k1 * CH4 * H2O
CO2 = CO * H2O / k2

eq = H2O + CO + 2*CO2 - 2*O
coeffs = get_coeffs(eq, CH4)
print("\nPolynomial coefficients for CH4 in HCO chemistry\n"
  "(Note that HL2016 makes an approximation for H2O (Eq (20)) that leads\n"
  "to different polynomial coefficients than these):")
for c in coeffs:
  print("  {},".format(c))
# Note that HL2016 makes an approximation for H2O (Eq (20)) that leads
# to different polynomial coefficients.


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# HCO chemistry with 6 molecules (no HCN, NH3, N2):

# Solve for CO:
CO  = s.symbols("CO")
H2O = k2 * (f*O-CO) / (k2+CO)
CO2 = CO * H2O / k2
CH4 = CO / (k1*H2O)
C2H2 = k3 * CH4**2
C2H4 = C2H2 / k4

# Polynomial coefficients sorted from lowest to highest degree:
eq = CH4 + CO + CO2 + 2*C2H2 + 2*C2H4 - f*C
coeffs = get_coeffs(eq, CO)[::-1]
print("\n\nPolynomial coefficients for CO in HCO chemistry with six species:")
for c in coeffs:
  print("  {},".format(c))


# Solve for H2O:
H2O = s.symbols("H2O")
CO  = k2 * (f*O-H2O) / (k2+H2O)
CO2 = CO * H2O / k2
CH4 = CO / (k1*H2O)
C2H2 = k3 * CH4**2
C2H4 = C2H2 / k4

# Polynomial coefficients sorted from lowest to highest degree:
eq = CH4 + CO + CO2 + 2*C2H2 + 2*C2H4 - f*C
coeffs = get_coeffs(eq, H2O)[::-1]
print("\n\nPolynomial coefficients for H2O in HCO chemistry with six species:")
for c in coeffs:
  print("  {},".format(c))


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# HCNO chemistry with 6-molecules (no CO2, C2H2, C2H4):

# Solve for CO:
CO  = s.symbols("CO")
H2O = 2*O - CO
CH4 = CO / (k1*H2O)
HCN = 2*C - CH4 - CO
NH3 = HCN / (CH4*k6)
N2  = k5 * NH3**2

# Polynomial coefficients sorted from lowest to highest degree:
eq = 2*N2 + NH3 + HCN - 2*N
coeffs = get_coeffs(eq, CO)[::-1]
print("\n\nPolynomial coefficients for CO in HCNO chemistry with six species:")
for c in coeffs:
  print("  {},".format(c))


# Solve for H2O:
H2O = s.symbols("H2O")
CO  = 2*O - H2O
CH4 = CO / (k1*H2O)
HCN = 2*C - CH4 - CO
NH3 = HCN / (CH4*k6)
N2  = k5 * NH3**2

# Polynomial coefficients sorted from lowest to highest degree:
eq = 2*N2 + NH3 + HCN - 2*N
coeffs = get_coeffs(eq, H2O)[::-1]
print("\n\nPolynomial coefficients for H2O in HCNO chemistry with six species:")
for c in coeffs:
  print("  {},".format(c))


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# HCNO chemistry 8-molecules (no CO2):

# Solve for CO:
CO = s.symbols("CO")
H2O  = f*O - CO
CH4  = CO / (k1*H2O)
C2H2 = k3 * CH4**2
C2H4 = C2H2 / k4
HCN  = f*C - CH4 - CO - 2*C2H2 - 2*C2H4
NH3  = HCN / (CH4*k6)
N2   = k5 * NH3**2

# Polynomial coefficients sorted from lowest to highest degree:
eq = 2*N2 + NH3 + HCN - f*N
coeffs = get_coeffs(eq, CO)[::-1]
scoeffs = s.factor(coeffs)
print("\n\nPolynomial coefficients for CO in HCNO chemistry with eight"
      " species:")
for c in scoeffs:
  print("  {},".format(c))

# Solve for H2O:
H2O = s.symbols("H2O")
CO  = f*O - H2O
CH4  = CO / (k1*H2O)
C2H2 = k3 * CH4**2
C2H4 = C2H2 / k4
HCN  = f*C - CH4 - CO - 2*C2H2 - 2*C2H4
NH3  = HCN / (CH4*k6)
N2   = k5 * NH3**2

# Polynomial coefficients sorted from lowest to highest degree:
eq = 2*N2 + NH3 + HCN - f*N
coeffs = get_coeffs(eq, H2O)[::-1]
scoeffs = s.factor(coeffs)
print("\n\nPolynomial coefficients for H2O in HCNO chemistry with eight"
      " species:")
for c in scoeffs:
  print("  {},".format(c))
