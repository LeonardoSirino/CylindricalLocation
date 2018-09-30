from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m

diameter = 400.0
height = 1000.0
semiperimeter = 240

Locate = CylindricalLocation(diameter, height)
Locate.setCalcMode('section')
Locate.set_semiPerimeter(semiperimeter)
Locate.SetVelocity(5)

Locate.StructuredSensorDistribution(
    lines=2, sensorsInLine=3, x0=0, y0=height / 3, dx=(diameter * m.pi) / 3, dy=height / 3, aligned=False)
Locate.StructuredSensorDistribution(
    lines=1, sensorsInLine=2, x0=0, y0=-semiperimeter / 2, dx=(diameter * m.pi) / 2, dy=0, aligned=False)
Locate.StructuredSensorDistribution(lines=1, sensorsInLine=2, x0=(
    diameter * m.pi) / 4, y0=height + semiperimeter / 2, dx=(diameter * m.pi) / 2, dy=0, aligned=False)


t = Locate.returnDeltaT(diameter * m.pi / 2, height / 2, [-1], True)

"""
data = []
i = 0
for AT in t:
    data.append((i, AT))
    i += 1

print(data)
"""
data = [(0, 96.67619767441275), (1, 20.198969522816547), (2, 20.19896952281656), (3, 56.8303869925099), (4, 0.0),
        (5, 56.8303869925099), (6, 137.3613916063625), (7, 90.6666833201405), (8, 104.2503305514595), (9, 104.25033055145903)]

print("Posição verdadeira:")
print("x: " + str(round(diameter * m.pi / 2, 4)) +
      " / y: " + str(round(height / 2, 4)))
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