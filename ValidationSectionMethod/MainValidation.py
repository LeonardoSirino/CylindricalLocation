from SectionClasses import *
import numpy as np

s = 50
xpos = []
s_sec = []
s_real = []
s_plan = []
for s in np.linspace(1, 50, num=20):
    for x in np.linspace(10, 314, num=100):

        xpos.append(x)

        point.f = 0.5
        point.diameter = 100
        point.divs = 100

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
        s_real.append(sreal / 2)

        calc = CalcSection()
        ssec = calc.distancePoints(point1, point2)
        s_sec.append(ssec)

        dplan1 = np.sqrt((x2 - x1)**2 + (s2 - s1)**2)
        dplan2 = np.sqrt((x2 - x1 + point.diameter * m.pi)**2 + (s2 - s1)**2)
        s_plan.append(min(dplan1, dplan2))

    s_plan.append(m.nan)
    s_real.append(m.nan)
    s_sec.append(m.nan)
    xpos.append(m.nan)

plt.plot(xpos, s_real, xpos, s_sec, xpos, s_plan)
plt.legend(['real', 'section', 'plan'])
plt.show()
