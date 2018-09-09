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
        self.mode = "reg"
        """Modos:
        reg: usa regressão para calcular arco
        inc: usa método incremental para calcular o arco
        """

    def RegPos(self, s):
        """ Regressão da função logística
        m = 0.84322062
        n = 1.43623984
        amp = 0.28950002

        R = (np.exp(s / m) - np.exp(n)) / (2 * amp * (np.exp(s / m) + np.exp(n)))
        """
        pol = [-0.04520616,  0.38323073, -1.36785798,  2.66208137, -3.10204898,  2.24771868,
               0.01275433, -1.00151578]
        R = np.polyval(pol, s)

        return R

    def AuxCoordsGeodesic(self):
        lon = self.x / (point.diameter * m.pi) * 360
        res = self.cap.Direct(lat1=0, lon1=lon, azi1=0, s12=self.s)
        self.lon = lon
        self.lat = res.get('lat2')

    def AuxCoordsSection(self):
        lon = self.x / (point.diameter * m.pi) * 360
        f = point.f
        a = point.diameter / 2

        if self.mode == "inc":
            sf = self.s
            s = 0
            R1 = a
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

        elif self.mode == "reg":
            R = self.RegPos(self.s / a)
            r = a * abs(R)
            try:
                z1 = a * f * m.sqrt(1 - r**2 / a**2)
            except ValueError:
                z1 = 0
                r = a

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
        """
        m = 0.82444066
        n = 1.46895599
        amp = 0.29505693
        s = m * (np.log((amp * R + 0.5) / (1 - amp * R - 0.5)) + n)
        """
        pol = [9.94406631e-01, -5.42331127e-13, -1.67276958e+00,  9.30333908e-13,
               9.92271096e-01, -4.52365676e-13, -1.60763495e-01,  8.69627507e-14,
               1.01106572e+00,  1.21105713e+00]
        s = np.polyval(pol, R)

        return s

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
