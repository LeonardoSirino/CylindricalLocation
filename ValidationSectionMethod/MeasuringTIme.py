from SectionClasses import *
import numpy as np
import time

tt_geo = 0
tt_section = 0
tt_plan = 0
s = 50
x = 100
N = 10000

point.f = 0.5
point.diameter = 100
point.divs = 100

t0 = time.time()
for i in range(0, N):
    x1 = x
    s1 = s
    point1 = point(x1, s1)
    point1.AuxCoordsGeodesic()

    x2 = 0
    s2 = s
    point2 = point(x2, s2)
    point2.AuxCoordsGeodesic()

    res = point1.cap.Inverse(lat1=point1.lat, lat2=point2.lat,
                                lon1=point1.lon, lon2=point2.lon)
    sreal = res.get("s12")

t1 = time.time()
tt_geo += t1-t0

t0 = time.time()
for i in range(0, N):

    calc = CalcSection()
    ssec = calc.distancePoints(point1, point2)

t1 = time.time()
tt_section += t1-t0

t0 = time.time()
for i in range(0, N):

    dplan1 = np.sqrt((x2 - x1)**2 + (s2 - s1)**2)
    dplan2 = np.sqrt((x2 - x1 + point.diameter * m.pi)**2 + (s2 - s1)**2)
    dplan = min(dplan1, dplan2)

t1 = time.time()
tt_plan += t1-t0

print("Tempo total geodésica:")
print(tt_geo)
print("Tempo total seção:")
print(tt_section)
print("Tempo total planificado:")
print(tt_plan)