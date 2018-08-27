import geographiclib.geodesic as geo
import numpy as np
import math as m
import time


class point():
    """Classe para agrupamento das propriedades de um ponto num elipsoide da revolução
    """
    diameter = 100
    f = 1 / 2
    divs = 100

    def __init__(self, x, s):
        self.x = x
        self.s = s
        self.a = point.diameter / 2
        self.cap = geo.Geodesic(point.diameter, point.f)

    def AuxCoordsGeodesic(self):
        lon = self.x / (point.diameter * m.pi) * 360
        res = self.cap.Direct(lat1=0, lon1=lon, azi1=0, s12=self.s)
        self.lon = lon
        self.lat = res.get('lat2')

    def AuxCoordsSection(self):
        sf = self.s
        s = 0
        a = point.diameter / 2
        R1 = a
        f = point.f
        z1 = 0
        dR = 2 * a / point.divs
        while s < sf:
            R2 = R1 - dR
            z2 = a * f * m.sqrt(1 - R2**2 / a**2)
            ds = m.sqrt((R2 - R1)**2 + (z2 - z1)**2)
            s += ds
            R1 = R2
            z1 = z2

        r = a - R1
        lon = self.x / (point.diameter * m.pi) * 360

        self.xcap = r * m.cos(m.radians(lon))
        self.ycap = r * m.sin(m.radians(lon))
        self.zcap = z1

    def __ellipseArc(self, Ri, Rf):
        s = 0
        a = point.diameter / 2
        R1 = a
        f = point.f
        z1 = a * f * m.sqrt(1 - Ri**2 / a**2)
        dR = (Rf - Ri) / point.divs
        radius = np.linspace(Ri, Rf, num=point.divs)
        for R in radius:
            z2 = a * f * m.sqrt(1 - R**2 / a**2)
            ds = m.sqrt(dR**2 + (z2 - z1)**2)
            s += ds
            z1 = z2


class CalcSection():
    """Métodos relacionados à distância por seccionamento
    """

    def __init__(self):
        pass

    def centerDistance(self, point1, point2):
        x1 = point1.xcap
        y1 = point1.ycap
        x2 = point2.xcap
        y2 = point2.ycap
        d = np.abs(x2 * y1 - y2 * x1) / np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        return d

    def reductionFactor(self, point1, point2):
        """Redução do diâmetro da elipse em função da sua distância do centro
        """
        a = point1.a
        d = self.centerDistance(point1, point2)
        redF = np.sqrt((a**2 - d**2) / a**2)

        return redF

    def distancePoints(self, point1, point2):
        point1.AuxCoordsSection()
        point2.AuxCoordsSection()
        redF = self.reductionFactor(point1, point2)
        a = redF * point1.a
        """Agora é preciso calcular a distância dos pontos para o centor dessa nova
        elipse, com isso se têm dados suficientes para calcular a distância entre os pontos
        """



point.f = 0.5
point.diameter = 100
point.divs = 100

x1 = 100
s1 = 20
point1 = point(x1, s1)
point1.AuxCoordsGeodesic()

x2 = 100
s2 = 60.5
point2 = point(x2, s2)
point2.AuxCoordsGeodesic()

res = point1.cap.Inverse(lat1=point1.lat, lat2=point2.lat,
                         lon1=point2.lon, lon2=point2.lon)
print("Distância verdadeira: ")
print(res.get("s12"))
