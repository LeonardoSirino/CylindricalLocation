from Routines.CilLoc import *

Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)
CoordList = [{"Xcord": 200.0, "Ycord": 1200.0},
             {"Xcord": 100.0, "Ycord": 900.0},
             {"Xcord": 200.0, "Ycord": 700.0},
             {"Xcord": 2500.0, "Ycord": -740 * m.pi / 4}]

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

Locate.calcAllDist(600, 900, True)