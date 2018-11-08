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

# Posição dos sensores
Locate.StructuredSensorDistribution(
    lines=2, sensorsInLine=3, x0=0, y0=height / 3, dx=(diameter * m.pi) / 3, dy=height / 3, aligned=False)
Locate.StructuredSensorDistribution(
    lines=1, sensorsInLine=2, x0=0, y0=-semiperimeter / 2, dx=(diameter * m.pi) / 2, dy=0, aligned=False)
Locate.StructuredSensorDistribution(lines=1, sensorsInLine=2, x0=(
    diameter * m.pi) / 4, y0=height + semiperimeter / 2, dx=(diameter * m.pi) / 2, dy=0, aligned=False)

# Ponto de teste
xp = diameter * m.pi * 0.7
yp = 990

t = Locate.returnDeltaT(xp, yp, [-1], 'geodesic')
data = []
i = 0
for AT in t:
    data.append((i, AT))
    i += 1

print("Posição verdadeira:")
print("x: " + str(round(xp, 4)) +
      " / y: " + str(round(yp, 4)))
print("\n")

t0 = time.time()

x = Locate.simpleLocation(data)

t1 = time.time()
print("Localização simplificada:")
print(x)
print("Tempo decorrido: " + str(t1 - t0))
print("\n")

t0 = time.time()

x = Locate.completeLocation(data)

t1 = time.time()
print("Localização completa:")
print(x)
print("Tempo decorrido: " + str(t1 - t0))

# Locate.fCostMap(data)