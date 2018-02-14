import geographiclib.geodesic as geo
import matplotlib.pyplot as plt


calc = geo.Geodesic(1, 1 / 2)
initLat = 0
initLon = 0
finalLat = 0
finalLon = 179
numDiv = 100

l = calc.InverseLine(lat1=initLat, lon1=initLon,
                     lat2=finalLat, lon2=finalLon)
increment = l.s13 / numDiv

lat = initLat
lon = initLon
s = 0
lats = []
lons = []
for i in range(numDiv+1):
    g = l.Position(s)
    lats.append(g.get("lat2"))
    lons.append(g.get("lon2"))
    s = min(s + increment, l.s13)

plt.plot(lons, lats)
plt.show()
