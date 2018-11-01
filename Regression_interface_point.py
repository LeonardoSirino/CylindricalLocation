from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import matplotlib.pyplot as plt
import time
import math as m

# Parâmetros do vaso
diameter = 400.0
f = 0.5


def InterfacePoint(alpha, rx_wall, rh_wall, rx_cap, rh_cap):

    # Configurações do algoritmo
    height = alpha * diameter
    Locate = CylindricalLocation(diameter, height)
    Locate.setCalcMode('section')
    Locate.set_f(f)
    semiperimeter = Locate.SemiPerimeter
    Locate.SetVelocity(5)

    # Ponto no tampo
    Pcap = VesselPoint(0, 0, 0)
    Pcap.Xcord = rx_cap * diameter * m.pi
    Pcap.Ycord = height + rh_cap * Locate.SemiPerimeter
    Pcap.OnCap = True
    Pcap.Cap = "sup"
    Locate._CylindricalLocation__AuxCoords(Pcap)

    # Ponto na parede
    Pwall = VesselPoint(0, 0, 1)
    Pwall.Xcord = rx_wall * diameter * m.pi
    Pwall.Ycord = height * rh_wall
    Pwall.OnCap = False

    dist, OtimP = Locate._CylindricalLocation__DistWalltoCap(Pcap, Pwall)
    x_trans = OtimP.Xcord / (diameter * m.pi)
    if x_trans > 1:
        x_trans -= 1
    elif x_trans < 0:
        x_trans += 1

    return x_trans

x = InterfacePoint(1.5, 0.1, 0.3, 0.5, 0.2)
print(x)