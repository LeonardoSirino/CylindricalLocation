import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def ArcReg(x, m, n):
    y = []
    amp = 0.40
    for pos in x:
        res = m * (np.log((amp * pos + 0.5)/(1 - amp*pos - 0.5)) + n)
        if res < 0:
            res = 0

        y.append(res)
    
    return y

def elipseArc(a, f):
    divs = 500
    s = 0
    Ri = -a
    Rf = a
    z1 = a * f * m.sqrt(1 - Ri**2 / a**2)
    dR = (Rf - Ri) / divs
    radius = np.linspace(Ri, Rf, num=divs)
    arc = []
    for R in radius:
        z2 = a * f * m.sqrt(1 - R**2 / a**2)
        ds = m.sqrt(dR**2 + (z2 - z1)**2)
        s += ds
        z1 = z2
        arc.append(s / a)

    pos = np.linspace(-1, 1, num=divs)

    fit, other = curve_fit(ArcReg, pos, arc, bounds=(0.01, [200., 2000.]))
    print(fit)
    y_reg = ArcReg(pos, fit[0], fit[1])
    residue = []
    for (real, adjust) in zip(arc, y_reg):
        residue.append(real - adjust)

    plt.plot(pos, residue)

    """
    plt.plot(pos, arc)
    plt.plot(pos, y_reg)
    plt.show()
    """

    return s, fit[0], fit[1]

s, v, n = elipseArc(1, 0.5)
"""
f_vec = np.linspace(0.001, 1, num=50)
v_vec = []
n_vec = []
for f in f_vec:
    s, v, n = elipseArc(1, f)
    v_vec.append(v)
    n_vec.append(n)

plt.show()
plt.plot(f_vec, v_vec, f_vec, n_vec)
plt.legend(["Parâmetro v", "Parâmetro n"])
plt.show()
"""