import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

f = 0.38

def elipseArc(a, f):
    divs = 100000
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

    fit = np.polyfit(pos, arc, order + 2)
    print("Regressão do arco em função da posição")
    print(fit)
    y_reg = np.polyval(fit, pos)

    residue = []
    for (real, adjust) in zip(arc, y_reg):
        residue.append(abs(real - adjust))

    print("Desvio máximo: " + str(np.max(residue)))

    # Vetores para plot
    index_plt = np.linspace(0, divs - 1, num = 11)
    pos_plt = []
    arc_plt = []
    for k in index_plt:
        k = int(k)
        pos_plt.append(pos[k])
        arc_plt.append(arc[k])

    index_plt = np.linspace(0, divs - 1, num = 10)
    posreg_plt = []
    yreg_plt = []
    for k in index_plt:
        k = int(k)
        posreg_plt.append(pos[k])
        yreg_plt.append(y_reg[k])

    plt.plot(pos_plt, arc_plt, '*')
    plt.plot(posreg_plt, yreg_plt, '.')
    plt.title("Regressão do arco")
    plt.ylabel("s / a")
    plt.xlabel("R / a")
    plt.legend(["Real", "Regressão"])
    plt.show()

    fit2 = np.polyfit(arc, pos, order)
    y_reg2 = np.polyval(fit2, arc)
    
    print("Regressão da posição em função do arco")
    print(fit2)
    residue2 = []
    for (real, adjust) in zip(pos, y_reg2):
        residue2.append(abs(real - adjust))

    print("Desvio máximo: " + str(np.max(residue2)))

    # Vetores para plot
    index_plt = np.linspace(0, divs - 1, num = 11)
    pos_plt = []
    arc_plt = []
    for k in index_plt:
        k = int(k)
        pos_plt.append(pos[k])
        arc_plt.append(arc[k])

    index_plt = np.linspace(0, divs - 1, num = 10)
    arcreg_plt = []
    yreg2_plt = []
    for k in index_plt:
        k = int(k)
        arcreg_plt.append(arc[k])
        yreg2_plt.append(y_reg2[k])

    plt.plot(arc_plt, pos_plt, '*')
    plt.plot(arcreg_plt, yreg2_plt, '.')
    plt.title("Regressão da posição")
    plt.xlabel("s / a")
    plt.ylabel("R / a")
    plt.legend(["Real", "Regressão"])
    plt.show()

    return s, fit[0], fit[1], fit[2]

order = 7
s, v, n, p = elipseArc(100, f)