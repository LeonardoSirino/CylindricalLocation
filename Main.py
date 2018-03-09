from Routines.CilLoc import *
import numpy as np
import time

diameter = 1000.0
height = 1000.0

Locate = CylindricalLocation(diameter, height, 750 * m.pi / 4)

""" Grid uniforme sobre o vaso de pressão"""
GridElements = 10
Xcords = np.linspace(0, diameter * m.pi, num=GridElements)
Ycords = np.linspace(1050, 1550, num=GridElements)

CoordList = []

for X in Xcords:
    for Y in Ycords:
        CoordList += [{"Xcord": X, "Ycord": Y}]


""" Pontos aleatórios
CoordList = [{"Xcord": 200.0, "Ycord": 1200.0},
             {"Xcord": 100.0, "Ycord": 900.0},
             {"Xcord": 200.0, "Ycord": 700.0},
             {"Xcord": 2500.0, "Ycord": -300}]
"""

""" Apenas um ponto
CoordList = [{"Xcord": 500.0, "Ycord": 1500.0}]
"""

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

t0 = time.time()

distances = Locate.calcAllDist(diameter * m.pi / 2, -500, True)

t1 = time.time()
print("Tempo decorrido: " + str(t1 - t0))

file = open("Output.txt", "w")
file.write("Tempo decorrido: " + str(t1 - t0) + "\n")
for dist in distances:
    file.write(str(dist) + "\n")
