from scipy.special import ellipkinc as elip
import numpy as np
import matplotlib.pyplot as plt
import math as m

divs = 10000
f = 0.5
ang = np.linspace(0, np.pi, num = divs)
y = []
for i in ang:
    res = elip(i, f)
    y.append(res)

s = 0
a = 1
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

plt.plot(ang, y, ang, arc)
plt.show()