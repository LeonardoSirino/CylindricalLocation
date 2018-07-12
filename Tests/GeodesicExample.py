import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import math as m

a = 1
calc = geo.Geodesic(a, 0.5)
initLat = 0
initLon = 0
finalLat = 0
finalLon = 175
numDiv = 100

l = calc.InverseLine(lat1=initLat, lon1=initLon,
                     lat2=finalLat, lon2=finalLon)
increment = l.s13 / numDiv

lat = initLat
lon = initLon
x_data = []
y_data = []
s = 0
lats = []
lons = []
for i in range(numDiv + 1):
    g = l.Position(s)
    lat = m.radians(g.get("lat2"))
    r = a * (1 - m.sin(lat))  # esse cálculo não está correto!!!
    lon = m.radians(g.get("lon2"))
    x = a - r * m.cos(lon)
    y = a - r * m.sin(lon)
    x_data.append(x)
    y_data.append(y)
    lats.append(g.get("lat2"))
    lons.append(g.get("lon2"))
    s = min(s + increment, l.s13)

plt.plot(x_data, y_data)
plt.show()
