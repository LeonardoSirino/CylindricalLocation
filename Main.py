from Routines.CilLoc import *
import numpy as np

Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)

CoordList = []
for height in np.linspace(-500, -10, num=5):
    CoordList += [{"Xcord": 200, "Ycord": height}]

for X in np.linspace(200, 3000, num=5):
    CoordList += [{"Xcord": X, "Ycord": -500}]

"""
CoordList = [{"Xcord": 200.0, "Ycord": 1200.0},
             {"Xcord": 100.0, "Ycord": 900.0},
             {"Xcord": 200.0, "Ycord": 700.0},
             {"Xcord": 2500.0, "Ycord": -300}]
"""

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

Locate.calcAllDist(1500, 1400, True)
