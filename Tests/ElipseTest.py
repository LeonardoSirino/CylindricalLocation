import math as m
import matplotlib.pyplot as plt
import numpy as np

N = 100


def calcElipseCoords(a, f, s0):
    s = 0
    x1 = a
    y1 = 0
    dx = 2 * a / N
    while s < s0:
        x2 = x1 - dx
        y2 = a * f * m.sqrt(1 - x2**2 / a**2)
        ds = m.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        s += ds
        x1 = x2
        y1 = y2

    return (x1, y1)


def ellipseArc(Ri, Rf, a, f):
    s = 0
    z1 = a * f * m.sqrt(1 - Ri**2 / a**2)
    dR = (Rf - Ri) / N
    radius = np.linspace(Ri, Rf, num = N)
    for R in radius:
        z2 = a * f * m.sqrt(1 - R**2 / a**2)
        ds = m.sqrt(dR**2 + (z2 - z1)**2)
        s += ds
        z1=z2

    return s

conv_s=[]
steps=range(10, 5000, 10)
for i in steps:
    N=i
    s=ellipseArc(0.1, 0.5, 1, 0.5)
    conv_s.append(s)


plt.plot(steps, conv_s)
plt.show()
