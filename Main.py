from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
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
Pcap = VesselPoint()
Pcap.Xcord = 0.2 * diameter * m.pi
Pcap.Ycord = height + 100
Pcap.OnCap = True
Pcao.Cap = "sup"

# Ponto na parede
Pwall = VesselPoint()
Pwall.Xcord = 0.3 * diameter * m.pi
Pwall.Ycord = height * 0.2
Pwall.OnCap = False


# Ponto de teste
xp = diameter * m.pi * 0.9
yp = 1100

t = Locate.returnDeltaT(xp, yp, [-1], 'geodesic')
print(t)