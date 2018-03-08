from Routines.CilLoc import *
import numpy as np
import time

Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)

GridElements = 7
Xcords = np.linspace(0,2800, num = GridElements)
Ycords = np.linspace(-500, 1500, num = GridElements)

CoordList= []

for X in Xcords:
    for Y in Ycords:
        CoordList += [{"Xcord": X, "Ycord": Y}]

""" Pontos aleatórios
CoordList = [{"Xcord": 200.0, "Ycord": 1200.0},
             {"Xcord": 100.0, "Ycord": 900.0},
             {"Xcord": 200.0, "Ycord": 700.0},
             {"Xcord": 2500.0, "Ycord": -300}]
"""

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

t0 = time.time()

distances = Locate.calcAllDist(2000,-100, True)

t1 = time.time()
print("Tempo decorrido: " + str(t1-t0))
