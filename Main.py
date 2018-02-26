from Routines.CilLoc import *

Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)
CoordList = [{"Xcord": 200.0, "Ycord": 1050.0},
             {"Xcord": 1000.0, "Ycord": 900.0},
             {"Xcord": 500.0, "Ycord": -740 * m.pi / 4}]

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

Locate.calcAllDist(200,950)