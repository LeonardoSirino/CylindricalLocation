from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import matplotlib.pyplot as plt
import time
import math as m

# Parâmetros do vaso
diameter = 400.0
height = 1000.0
f = 0.5

# Configurações do algoritmo
Locate = CylindricalLocation(diameter, height)
Locate.setCalcMode('section')
Locate.set_f(f)
semiperimeter = Locate.SemiPerimeter
Locate.SetVelocity(5)

# Ponto no tampo
Pcap = VesselPoint(0, 0, 0)
Pcap.Xcord = 0.5 * diameter * m.pi
Pcap.Ycord = height + 100
Pcap.OnCap = True
Pcap.Cap = "sup"
Locate._CylindricalLocation__AuxCoords(Pcap)

# Ponto na parede
Pwall = VesselPoint(0, 0, 1)
Pwall.Xcord = 0.3 * diameter * m.pi
Pwall.Ycord = height * 0.2
Pwall.OnCap = False

num_Xwall = 200
num_Ywall = 10
X_Pwall = np.linspace(0, diameter * m.pi, num=num_Xwall)
Y_Pwall = np.linspace(0, height, num=num_Ywall)
x_otim = np.zeros(num_Xwall)
legend = []

j = 0
for y in Y_Pwall:
    i = 0
    for x in X_Pwall:
        Pwall.Xcord = x
        Pwall.Ycord = y
        dist, OtimP = Locate._CylindricalLocation__DistWalltoCap(Pcap, Pwall)
        x_trans = OtimP.Xcord / (diameter * m.pi)
        if x_trans > 1:
            x_trans -= 1
        elif x_trans < 0:
            x_trans += 1
        x_otim[i] = x_trans

        i += 1

    plt.plot(X_Pwall / (diameter * m.pi), x_otim)
    legend.append(str(y / height))
    j += 1

plt.legend(legend)
plt.show()
