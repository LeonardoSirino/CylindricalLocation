import geographiclib.geodesic as geo
import numpy as np
import math as m
import time
import matplotlib.pyplot as plt


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
        self.cap = geo.Geodesic(point.diameter / 2, point.f)

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

        r = R1
        lon = self.x / (point.diameter * m.pi) * 360

        self.xcap = r * m.cos(m.radians(lon))
        self.ycap = r * m.sin(m.radians(lon))
        self.zcap = z1


class CalcSection():
    """Métodos relacionados à distância por seccionamento
    """

    def __init__(self):
        self.mode = "reg"
        """Modos:
        reg: usa regressão para calcular arco
        inc: usa método incremental para calcular o arco
        """

    def RegArc(self, R):
        m = 0.81010676
        n = 1.49757349
        amp = 0.3
        res = m * (np.log((amp * R + 0.5)/(1 - amp * R - 0.5)) + n)

        return res

    def RegPos(self, s):
        R = 0
        """Função para retornar a posição de um ponto do qual se sabe o comprimento do arco
        Ajustar uma função logística para essa curva
        """
        return R

    def centerLineDistance(self, point1, point2):
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
        d = self.centerLineDistance(point1, point2)
        redF = np.sqrt((a**2 - d**2) / a**2)

        return redF, d

    def elipseArc(self, a, u1, u2):
        if self.mode == "inc":
            f = point.f
            s = 0
            Ri = min([u1, u2])
            Rf = max([u1, u2])
            z1 = a * f * m.sqrt(1 - Ri**2 / a**2)
            dR = (Rf - Ri) / point.divs
            radius = np.linspace(Ri, Rf, num=point.divs)
            for R in radius:
                z2 = a * f * m.sqrt(1 - R**2 / a**2)
                ds = m.sqrt(dR**2 + (z2 - z1)**2)
                s += ds
                z1 = z2
        
        elif self.mode == "reg":
            Ri = min([u1, u2])
            Rf = max([u1, u2])
            si = self.RegArc(Ri / a)
            sf = self.RegArc(Rf / a)
            s = (sf - si) * a

        return s

    def distancePoints(self, point1, point2):
        point1.AuxCoordsSection()
        point2.AuxCoordsSection()
        redF, d = self.reductionFactor(point1, point2)
        r1q = point1.xcap**2 + point1.ycap**2
        u1 = m.sqrt(r1q - d**2)
        r2q = point2.xcap**2 + point2.ycap**2
        u2 = m.sqrt(r2q - d**2)

        v1c = np.array([-point1.xcap, -point1.ycap])
        v2c = np.array([-point2.xcap, -point2.ycap])
        v12 = np.array([point2.xcap - point1.xcap, point2.ycap - point1.ycap])

        theta1 = np.arccos(np.dot(v1c, v12) /
                           (np.linalg.norm(v1c) * np.linalg.norm(v12)))

        theta2 = np.arccos(np.dot(v2c, v12) /
                           (np.linalg.norm(v2c) * np.linalg.norm(v12)))

        if theta1 > m.pi / 2:
            u1 = -u1

        if theta2 > m.pi / 2:
            u2 = -u2

        a = redF * point1.a
        s = self.elipseArc(a, u1, u2)
        return s