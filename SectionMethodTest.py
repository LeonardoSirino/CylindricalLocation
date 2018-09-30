from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m

diameter = 400.0
height = 1000.0
semiperimeter = 240

Locate = CylindricalLocation(diameter, height)
Locate.setCalcMode("section")
Locate.set_semiPerimeter(semiperimeter)
Locate.SetVelocity(5)

Locate.AddSensor(500, -100)