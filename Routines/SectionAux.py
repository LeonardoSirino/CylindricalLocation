import math as m
import numpy as np
from scipy.optimize import curve_fit

def elipseArcReg(f, order, divs):
    s = 0
    Ri = -1
    Rf = 1
    z1 = f * m.sqrt(1 - Ri**2)
    dR = (Rf - Ri) / divs
    radius = np.linspace(Ri, Rf, num=divs)
    arc = []
    for R in radius:
        z2 = f * m.sqrt(1 - R**2)
        ds = m.sqrt(dR**2 + (z2 - z1)**2)
        s += ds
        z1 = z2
        arc.append(s)

    pos = np.linspace(-1, 1, num=divs)

    fit = np.polyfit(pos, arc, order + 2)
    fit2 = np.polyfit(arc, pos, order)

    return fit, fit2

"""
order = 7
f = 0.5
divs = 100000
f1, f2 = elipseArcReg(f, order, divs)
print(f1)
print(f2)
"""