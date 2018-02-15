import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as opt
import numpy as np

diameter = 100
f = 1 / 2
semiPerim = m.pi * (diameter * (1 + f) / 2) / 4
print("Semi-Perimetro: " + str(semiPerim))
# Point 1
P1 = {"Xcord": 10, "Ycord": 10}
# Point 2
P2 = {"Xcord": 20, "Ycord": -40}

tampo = geo.Geodesic(diameter, f)

# Determining latitude and longitude for Point 2
P2lon = P2.get("Xcord") / diameter * 360 - 180
print("P2lon: " + str(P2lon))
P2Geo = tampo.Direct(lat1=0, lon1=P2lon, s12=abs(P2.get("Ycord")), azi1=0)
P2lat = P2Geo.get("lat2")
print("P2lat: " + str(P2lat))


def CalCDist(diamPos, P1, P2lat, P2lon, semiP, diam):
    d1 = m.sqrt((diamPos - P1.get("Xcord"))**2 + P1.get("Ycord")**2)
    #print("d1: " + str(d1))
    AuxLon = diamPos / diam * 360 - 180
    #print("auxlon: " + str(AuxLon))
    distTampo = tampo.Inverse(lat1=0, lat2=P2lat, lon1=AuxLon, lon2=P2lon)
    #print("distTampo: " + str(distTampo.get("s12")))
    return d1 + distTampo.get("s12")


def Fopt(auxdiam):
    return CalCDist(auxdiam, P1, P2lat, P2lon, semiPerim, diameter)


res = opt.minimize(Fopt, 10, method="nelder-mead",
                   options={"xtol": 1e-8, "disp": True})
print(res.x)

AuxPoints = np.linspace(P1.get("Xcord"), P2.get("Xcord"), 50)
Distances = []
for AuxPoint in AuxPoints:
    dist = CalCDist(AuxPoint, P1, P2lat, P2lon, semiPerim, diameter)
    Distances.append(dist)

plt.plot(AuxPoints, Distances)
plt.show()
