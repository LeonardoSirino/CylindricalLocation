import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as opt
import numpy as np

diameter = 100
f = 1 / 2
height = 100
semiPerim = m.pi * (diameter * (1 + f) / 2) / 4
print("Semi-Perimetro: " + str(semiPerim))
# Point 1
P1 = {"Xcord": 10, "Ycord": 50}
# Point 2
P2 = {"Xcord": 250, "Ycord": -40}

tampo = geo.Geodesic(diameter, f)

# Determining latitude and longitude for Point 2
P2lon = P2.get("Xcord") / (diameter * m.pi) * 360 - 180
print("P2lon: " + str(P2lon))
P2Geo = tampo.Direct(lat1=0, lon1=P2lon, s12=abs(P2.get("Ycord")), azi1=0)
P2lat = P2Geo.get("lat2")
print("P2lat: " + str(P2lat))


def CalCDist(diamPos, P1, P2lat, P2lon, semiP, diam):
    d1 = m.sqrt((diamPos - P1.get("Xcord"))**2 + P1.get("Ycord")**2)
    #print("d1: " + str(d1))
    AuxLon = diamPos / diam * m.pi * 360 - 180
    #print("auxlon: " + str(AuxLon))
    distTampo = tampo.Inverse(lat1=0, lat2=P2lat, lon1=AuxLon, lon2=P2lon)
    #print("distTampo: " + str(distTampo.get("s12")))
    return d1 + distTampo.get("s12")


def Fopt(auxdiam):
    return CalCDist(auxdiam, P1, P2lat, P2lon, semiPerim, diameter)


res = opt.minimize_scalar(Fopt, bounds=(
    P1.get("Xcord"), P2.get("Xcord")), method="bounded")
print("Solução, posição em x: " + str(res.x))
lonMin = res.x / (diameter * m.pi) * 360 - 180
print("Solução, longitude: " + str(lonMin))

Geoline = tampo.InverseLine(lat1=0, lat2=P2lat, lon1=lonMin, lon2=P2lon)
smax = tampo.Inverse(lat1=0, lat2=P2lat, lon1=lonMin, lon2=P2lon).get("s12")
print("smax: " + str(smax))

ycords = [P1.get("Ycord"), 0]
xcords = [P1.get("Xcord"), res.x]
for s in np.linspace(start=0, stop=smax):
    pos = Geoline.Position(s)
    posDiam = (pos.get("lon2") + 180) * diameter * m.pi / 360
    lataux = pos.get("lat2")
    Saux = tampo.Inverse(lat1=0, lon1=0, lat2=lataux, lon2=0).get("s12")
    ycords.append(-Saux)
    xcords.append(posDiam)

maxX = diameter * m.pi
VesselX = [0, maxX, maxX, 0, 0, 0, maxX, maxX]
VesselY = [-semiPerim, -semiPerim, 0, 0, -semiPerim, height, height, 0]


plt.plot(xcords, ycords, "g-")
plt.plot(VesselX, VesselY)
plt.plot([P1.get("Xcord"), P2.get("Xcord")], [P1.get("Ycord"), P2.get("Ycord")], "r")
plt.ylabel("Altura do vaso")
plt.xlabel("Distância a partir da geratriz")
plt.title("Caminho direto x caminho real")
plt.show()
plt.clf()


AuxPoints = np.linspace(P1.get("Xcord"), P2.get("Xcord"), 100)
Distances = []
for AuxPoint in AuxPoints:
    dist = CalCDist(AuxPoint, P1, P2lat, P2lon, semiPerim, diameter)
    Distances.append(dist)

plt.plot(AuxPoints, Distances)
plt.show()
