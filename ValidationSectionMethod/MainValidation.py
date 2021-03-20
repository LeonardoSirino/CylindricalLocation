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
arcs = np.linspace(0, 60.5, num=10)
max_errors = []
max_errors_plan = []
diameter = 100
for s in arcs:
    max_error = 0
    max_error_plan = 0
    for x in np.linspace(1, 314, num=500):

        xpos.append(x / diameter)

        point.f = 0.5
        point.diameter = diameter
        point.divs = 500

        t0 = time.time()

        x1 = x
        s1 = s
        point1 = point(x1, s1)
        point1.AuxCoordsGeodesic()

        x2 = 0
        s2 = 0
        point2 = point(x2, s2)
        point2.AuxCoordsGeodesic()

        res = point1.cap.Inverse(lat1=point1.lat, lat2=point2.lat,
                                 lon1=point1.lon, lon2=point2.lon)
        sreal = res.get("s12")

        t1 = time.time()
        tt_geo += t1 - t0

        s_real.append(sreal / diameter)

        t0 = time.time()

        calc = CalcSection()
        ssec = calc.distancePoints(point1, point2)

        t1 = time.time()
        tt_section += t1 - t0

        s_sec.append(ssec / diameter)

        t0 = time.time()

        dplan1 = np.sqrt((x2 - x1)**2 + (s2 - s1)**2)
        dplan2 = np.sqrt((x2 - x1 + point.diameter * m.pi)**2 + (s2 - s1)**2)
        dplan = min(dplan1, dplan2)

        t1 = time.time()
        tt_plan += t1 - t0

        s_plan.append(dplan / diameter)

        erro_plan.append(abs(sreal - dplan) / diameter)
        erro_section.append(abs(sreal - ssec) / diameter)

        if abs(sreal - ssec) / diameter > abs(max_error):
            max_error = abs((sreal - ssec) / diameter)

        if abs(dplan - sreal) / diameter > abs(max_error_plan):
            max_error_plan = abs((dplan - sreal) / diameter)

    max_errors.append(max_error)
    max_errors_plan.append(max_error_plan)
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
plt.legend(['Real', 'Seccionamento', 'Planificado'])
plt.title("Distâncias")
plt.xlabel('Posição X normalizada')
plt.ylabel('Distância normalizada')
plt.show()

plt.plot(xpos, erro_plan, xpos, erro_section)
plt.legend(['Erro planificado', 'Erro seccionamento'])
plt.title("Erros")
plt.xlabel('Posição X normalizada')
plt.ylabel("Erro normalizado")
plt.show('Erro normalizado')

# plt.plot(arcs, max_errors, arcs, max_errors_plan)
plt.plot(arcs, max_errors, arcs, max_errors_plan)
plt.legend(['Erro seccionamento', 'Erro planificado'])
plt.title("Erro máximo normalizado")
plt.xlabel("Distância normalizada até o corpo")
plt.ylabel("Erro normalizado")
plt.show()