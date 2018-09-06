from SectionClasses import *
import numpy as np
import time

xpos = []
s_sec = []
s_real = []
s_plan = []
erro_plan = []
erro_section = []
tt_geo = 0
tt_section = 0
tt_plan = 0
for s in np.linspace(10, 50, num=5):
    for x in np.linspace(1, 314, num=100):

        xpos.append(x)

        point.f = 0.5
        point.diameter = 100
        point.divs = 100

        t0 = time.time()

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

        s_real.append(sreal)

        t0 = time.time()

        calc = CalcSection()
        ssec = calc.distancePoints(point1, point2)

        t1 = time.time()
        tt_section += t1-t0

        s_sec.append(ssec)

        t0 = time.time()

        dplan1 = np.sqrt((x2 - x1)**2 + (s2 - s1)**2)
        dplan2 = np.sqrt((x2 - x1 + point.diameter * m.pi)**2 + (s2 - s1)**2)
        dplan = min(dplan1, dplan2)

        t1 = time.time()
        tt_plan += t1-t0

        s_plan.append(dplan)

        erro_plan.append(sreal - dplan)
        erro_section.append(sreal - ssec)

    s_plan.append(m.nan)
    s_real.append(m.nan)
    s_sec.append(m.nan)
    erro_plan.append(m.nan)
    erro_section.append(m.nan)
    xpos.append(m.nan)

print("Tempo total geodésica:")
print(tt_geo)
print("Tempo total seção:")
print(tt_section)
print("Tempo total planificado:")
print(tt_plan)
plt.plot(xpos, s_real, xpos, s_sec, xpos, s_plan)
plt.legend(['real', 'section', 'plan'])
plt.title("Distâncias")
plt.show()
plt.plot(xpos, erro_plan, xpos, erro_section)
plt.legend(['erro planificado', 'erro section'])
plt.title("Erros")
plt.show()