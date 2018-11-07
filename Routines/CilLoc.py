import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math as m
import numpy as np
import copy
import time

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

class VesselPoint():
    """Classe para definir as propriedades de um ponto no vaso
    """

    def __init__(self, Xcord, Ycord, ID):
        self.Xcord = Xcord
        self.Ycord = Ycord
        self.ID = ID

    def SetXcord(self, Xcord):
        self.Xcord = Xcord

    def SetLon(self, Lon):
        self.Lon = Lon

    def SetLat(self, Lat):
        self.Lat = Lat

    def SetOnCap(self, OnCap):
        self.OnCap = OnCap

    def SetCap(self, Cap):
        self.Cap = Cap

    def SetXcap(self, X):
        self.Xcap = X

    def SetYcap(self, Y):
        self.Ycap = Y

    def SetZcap(self, Z):
        self.Zcap = Z

    def __str__(self):
        text = "ID: " + \
            str(self.ID) + " x: " + str(self.Xcord) + " y: " + str(self.Ycord)
        return text


class CylindricalLocation():
    def __init__(self, diameter, height):
        self.diameter = diameter
        self.height = height
        self.SensorList = []
        self.veloc = 1
        self.__SensorID = -1
        self.__tempSensorList = None
        self.CalcMode = 'geodesic'
        """Modos:
        geodesic - usando biblioteca do Python - GeoplotLib
        section - usando seccionamento do tampo
        """
        self.SectionMode = 'reg'
        """Modos:
        reg: usa regressão para calcular arco
        inc: usa método incremental para calcular o arco
        """
        self.__ellipseDivs = 500
        self.__DivsTolerance = 100

        # Inicialização dos tempos acumulados
        self.t_samecap = 0
        self.t_wallcap = 0
        self.t_wall = 0
        self.t_captocap = 0
        self.i_samecap = 0
        self.i_wallcap = 0
        self.i_wall = 0
        self.i_captocap = 0

    def setCalcMode(self, mode):
        self.CalcMode = mode
        """Modos:
        geodesic - usando biblioteca do Python - GeoplotLib
        section - usando seccionamento do tampo
        """

    def setSectionMode(self, mode):
        self.SectionMode = mode
        """Modos:
        reg: usa regressão para calcular arco
        inc: usa método incremental para calcular o arco
        """

    def set_f(self, f):
        self.f = f
        self.cap = geo.Geodesic(self.diameter / 2, f)
        result = self.cap.Inverse(lat1=0, lon1=0, lat2=90, lon2=0)
        self.SemiPerimeter = result.get("s12")

    def set_semiPerimeter(self, SemiPerimeter):
        """Definição do valor do semiperímetro e calculo do valor de f correspondente

        Arguments:
            SemiPerimeter {[float]} -- [medida do semiperímetro]

        Returns:
            f[float] -- [razão de achatamento]
        """

        self.SemiPerimeter = SemiPerimeter

        def CalcSemiPerimeter(f):
            cap = geo.Geodesic(self.diameter / 2, f)
            result = cap.Inverse(lat1=0, lon1=0, lat2=90, lon2=0)
            CalcSP = result.get('s12')
            return CalcSP

        res = opt.minimize(lambda x: (CalcSemiPerimeter(
            x) - SemiPerimeter)**2, bounds=[(0, 0.999)], method='L-BFGS-B', x0=0.5)
        self.f = res.get("x")[0]
        self.cap = geo.Geodesic(self.diameter / 2, self.f)

    def SetVelocity(self, velocity):
        "Definição da velocidade em mm/s"
        self.veloc = velocity

    def PrintAllSensors(self):
        print("Sensores originais:")
        for sensor in self.SensorList:
            print(sensor)

    def __getSensorbyID(self, ID):
        for sensor in self.SensorList:
            if sensor.ID == ID:
                return sensor

        return None

    def __sectionPos(self, s):
        if self.SectionMode == "reg":
            a = self.diameter / 2
            f = self.f
            s = s / a
            pol = [-0.04520616,  0.38323073, -1.36785798,  2.66208137, -3.10204898,  2.24771868,
                   0.01275433, -1.00151578]
            R = np.polyval(pol, s) * a
            try:
                z = a * f * m.sqrt(1 - R**2 / a**2)
            except ValueError:
                if abs(a - abs(R)) < (a / self.__DivsTolerance):
                    R = a
                    z = 0
                else:
                    print("SectionPos")
                    print("a: " + str(a))
                    print("s: " + str(s))
                    print("R: " + str(R))
                    z = np.nan
                    R = np.nan
        elif self.SectionMode == "inc":
            sf = s
            f = self.f
            a = self.diameter / 2
            R1 = a
            s = 0
            z1 = 0
            dR = 2 * a / self.__ellipseDivs
            while s < sf:
                R2 = R1 - dR
                z2 = a * f * m.sqrt(1 - R2**2 / a**2)
                ds = m.sqrt((R2 - R1)**2 + (z2 - z1)**2)
                s += ds
                R1 = R2
                z1 = z2

            R = R1
            z = z1
        else:
            print("Modo inexistente")

        return (R, z)

    def __sectionArc(self, a, R1, R2):
        Ri = min([R1, R2])
        Rf = max([R1, R2])
        f = self.f

        if self.SectionMode == "reg":
            Ri = Ri / a
            Rf = Rf / a
            pol = [9.94406631e-01, -5.42331127e-13, -1.67276958e+00,  9.30333908e-13,
                   9.92271096e-01, -4.52365676e-13, -1.60763495e-01,  8.69627507e-14,
                   1.01106572e+00,  1.21105713e+00]
            si = np.polyval(pol, Ri) * a
            sf = np.polyval(pol, Rf) * a
            s = sf - si
        elif self.SectionMode == "inc":
            s = 0
            z1 = a * f * m.sqrt(1 - Ri**2 / a**2)
            dR = (Rf - Ri) / self.__ellipseDivs
            radius = np.linspace(Ri, Rf, num=self.__ellipseDivs)
            for R in radius:
                z2 = a * f * m.sqrt(1 - R**2 / a**2)
                ds = m.sqrt(dR**2 + (z2 - z1)**2)
                s += ds
                z1 = z2
        else:
            print("Modo inexistente")

        return s

    def __centerLineDistance(self, point1, point2):
        a = self.diameter / 2
        N = self.__DivsTolerance
        x1 = point1.Xcap
        y1 = point1.Ycap
        x2 = point2.Xcap
        y2 = point2.Ycap
        if abs(x1 - x2) < a / N and abs(y1 - y2) < a / N:
            # Evita erros no cálculo da distância até o centro quando P1 e P2 estão muito próximos
            d = a
        else:
            d = np.abs(x2 * y1 - y2 * x1) / \
                np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        return d

    def __reductionFactor(self, point1, point2):
        """Redução do diâmetro da elipse em função da sua distância do centro
        """
        a = self.diameter / 2
        d = self.__centerLineDistance(point1, point2)
        redF = np.sqrt((a**2 - d**2) / a**2)

        if np.isnan(redF) or a == d:
            redF = 0

        return redF, d

    def __AuxCoords(self, Coords):
        """
        Função para calcular as coordenadas auxliares (latitude e longitude) quando a coordenada estiver no tampo
        """
        if Coords.Ycord >= self.height or Coords.Ycord <= 0:
            if self.CalcMode == 'geodesic':
                """Calculo de variáveis auxiliares necessárias para o uso da geodésicas
                """

                lon = Coords.Xcord / (self.diameter * m.pi) * 360 - 180
                if Coords.Ycord >= self.height:
                    Coords.SetCap("sup")
                    s12 = Coords.Ycord - self.height
                else:
                    s12 = abs(Coords.Ycord)
                    Coords.SetCap("inf")

                EndCap = self.cap.Direct(lat1=0, lon1=lon, s12=s12, azi1=0)
                lat = EndCap.get("lat2")
                Coords.SetLat(lat)
                Coords.SetLon(lon)
            elif self.CalcMode == 'section':
                """Cálculo de variáveis auxiliares para o uso do seccionamento do tampo
                """
                lon = Coords.Xcord / (self.diameter * m.pi) * 2 * m.pi - m.pi
                if Coords.Ycord >= self.height:
                    Coords.SetCap("sup")
                    s = Coords.Ycord - self.height
                else:
                    s = abs(Coords.Ycord)
                    Coords.SetCap("inf")

                (R, z) = self.__sectionPos(s)
                R = abs(R)

                Coords.Xcap = R * m.cos(lon)
                Coords.Ycap = R * m.sin(lon)
                Coords.Zcap = z
                """
                print("Coordenada auxiliar tampo X: " + str(Coords.Xcap))
                print("Coordenada auxiliar tampo Y: " + str(Coords.Ycap))
                print("Coordenada auxiliar tampo Z: " + str(Coords.Zcap))
                """
            else:
                print("Modo inválido")

            Coords.SetOnCap(True)
        else:
            Coords.SetOnCap(False)

    def __AllSensorsAuxCoords(self):
        for sensor in self.SensorList:
            self.__AuxCoords(sensor)

    def __calcDist(self, P1, P2):
        """
        Função para calcular distância entre dois pontos em posições quaisquer do vaso
        """
        if P1.OnCap and P2.OnCap:
            if P1.Cap == P2.Cap:  # Dois pontos no mesmo tampo
                t0 = time.time()
                dist = self.__DistSameCap(P1, P2)
                t1 = time.time()
                self.t_samecap += t1 - t0
                self.i_samecap += 1
            else:  # Pontos em tampos opostos
                if P1.Cap == "sup":
                    Psup = P1
                    Pinf = P2
                else:
                    Psup = P2
                    Pinf = P1
                t0 = time.time()
                dist = self.__DistCaptoCap(Psup, Pinf)
                t1 = time.time()
                self.t_captocap += t1 - t0
                self.i_captocap += 1
        # Distância entre um ponto no casco e outro no tampo
        elif P1.OnCap ^ P2.OnCap:
            t0 = time.time()
            dist = self.__DistWalltoCap(P1, P2)
            t1 = time.time()
            self.t_wallcap += t1 - t0
            self.i_wallcap += 1
        else:  # Distância entre pontos no casco
            t0 = time.time()
            dist1 = np.sqrt((P1.Xcord - P2.Xcord) **
                           2 + (P1.Ycord - P2.Ycord)**2)
            # Clone à direita
            dist2 = np.sqrt((P1.Xcord - P2.Xcord + self.diameter *
                            m.pi) ** 2 + (P1.Ycord - P2.Ycord)**2)
            # Clone à esquerda
            dist3 = np.sqrt((P1.Xcord - P2.Xcord - self.diameter *
                            m.pi) ** 2 + (P1.Ycord - P2.Ycord)**2)

            dist = np.min([dist1, dist2, dist3])

            t1 = time.time()
            self.t_wall += t1 - t0
            self.i_wall += 1

        return dist

    def __DistSameCap(self, P1, P2):
        if self.CalcMode == 'geodesic':
            res = self.cap.Inverse(
                lat1=P1.Lat, lat2=P2.Lat, lon1=P1.Lon, lon2=P2.Lon)
            dist = res.get("s12")

        elif self.CalcMode == 'section':
            redF, d = self.__reductionFactor(P1, P2)
            if redF != 0:
                r1q = P1.Xcap**2 + P1.Ycap**2
                try:
                    u1 = m.sqrt(r1q - d**2)
                except ValueError:
                    if abs(r1q - d) < self.diameter / self.__DivsTolerance:
                        u1 = 0
                    else:
                        u1 = np.nan

                r2q = P2.Xcap**2 + P2.Ycap**2
                try:
                    u2 = m.sqrt(r2q - d**2)
                except ValueError:
                    if abs(r2q - d) < self.diameter / self.__DivsTolerance:
                        u2 = 0
                    else:
                        u2 = np.nan

                v1c = np.array([-P1.Xcap, -P1.Ycap])
                v2c = np.array([-P2.Xcap, -P2.Ycap])
                v12 = np.array([P2.Xcap - P1.Xcap, P2.Ycap - P1.Ycap])

                cosTheta1 = np.dot(v1c, v12) / \
                    (np.linalg.norm(v1c) * np.linalg.norm(v12))

                cosTheta2 = np.dot(v2c, v12) / \
                    (np.linalg.norm(v2c) * np.linalg.norm(v12))

                if cosTheta1 <= 0:
                    # Ponto está do outro lado do centro da elipse
                    u1 = -u1

                if cosTheta2 <= 0:
                    # Ponto está do outro lado do centro da elipse
                    u2 = -u2

                a = redF * self.diameter / 2
                dist = self.__sectionArc(a, u1, u2)
            else:
                dist = 0
        else:
            print("Modo inexistente")

        return dist

    def __DistWalltoCap(self, P1, P2):
        if P1.OnCap:
            Pcap = P1
            Pwall = P2
        else:
            Pcap = P2
            Pwall = P1

        if Pcap.Cap == "sup":
            Yaux = self.height
        else:
            Yaux = 0

        AuxPoint = VesselPoint(0, Yaux, -2)

        def CalcWallToCap(Xaux):
            """[Função para cálculo de distância entre um ponto no casco e outro no tampo]

            Arguments:
                Xaux {[float]} -- [Posição x do ponto auxiliar: ponto de transição casco-tampo]

            Returns:
                [totalDist] -- [Distância total entre pontos]
            """

            AuxPoint.SetXcord(Xaux)
            self.__AuxCoords(AuxPoint)
            AuxPoint.SetOnCap(False)

            t0 = time.time()

            dist1 = self.__calcDist(Pwall, AuxPoint)

            t1 = time.time()
            dt1 = t1 - t0
            self.t_wall -= dt1
            self.i_wall -= 1

            AuxPoint.SetOnCap(True)

            t0 = time.time()

            dist2 = self.__calcDist(AuxPoint, Pcap)

            t1 = time.time()
            dt2 = t1 - t0
            self.t_samecap -= dt2
            self.i_samecap -= 1

            totalDist = dist1 + dist2
            return totalDist

        InitGuess = opt.brute(CalcWallToCap, ((0, m.pi * self.diameter),), Ns=6)
        FinalSearch = opt.minimize(CalcWallToCap, x0=InitGuess, method="BFGS", options={'maxiter': 5})
        # print(FinalSearch) -- Resultado da minimização
        dist = FinalSearch.get("fun")

        return dist

    def __DistCaptoCap(self, Psup, Pinf):
        AuxPoint1 = VesselPoint(0, self.height, -2)
        AuxPoint2 = VesselPoint(0, 0, -2)

        def CalcCaptoCap(Xaux):
            """[Função para calcular a distância entre pontos em tampos opostos]

            Arguments:
                Xaux {[float array]} -- [vetor com as coordenadas x dos pontos auxiliares]

            Returns:
                [dist] -- [distância entre os pontos]
            """

            AuxPoint1.SetXcord(Xaux[0])
            AuxPoint2.SetXcord(Xaux[1])
            self.__AuxCoords(AuxPoint1)
            self.__AuxCoords(AuxPoint2)
            AuxPoint1.SetOnCap(True)
            dist1 = self.__calcDist(Psup, AuxPoint1)
            AuxPoint1.SetOnCap(False)
            AuxPoint2.SetOnCap(False)
            dist2 = self.__calcDist(AuxPoint1, AuxPoint2)
            AuxPoint2.SetOnCap(True)
            dist3 = self.__calcDist(AuxPoint2, Pinf)
            totalDist = dist1 + dist2 + dist3
            return totalDist

        # Chute inicial com otimização bruta
        SearchRange = (0, m.pi * self.diameter)
        InitGuess = opt.brute(CalcCaptoCap, (SearchRange, SearchRange), Ns=6)
        FinalSearch = opt.minimize(CalcCaptoCap, x0=InitGuess, method="BFGS")
        # print(FinalSearch)  # -- Resultado da minimização
        dist = FinalSearch.get("fun")
        MinPosSup = FinalSearch.get("x")[0]
        MinPosInf = FinalSearch.get("x")[1]
        AuxPoint1.SetXcord(MinPosSup)
        AuxPoint2.SetXcord(MinPosInf)

        return dist

    def __DistVClone(self, Source, Sensor):
        # Mudar o calculo para uma condição simplificada, assumindo que o caminho da onda é vertical até a interface, dessa forma não haverá minimização envolvida
        SemiHeight = self.height / 2
        dx1 = abs(Source.Xcord - Sensor.Xcord)
        dx2 = abs(Source.Xcord - Sensor.Xcord + m.pi * self.diameter)
        dx3 = abs(Source.Xcord - Sensor.Xcord - m.pi * self.diameter)
        dx = np.min([dx1, dx2, dx3])
        dy = Sensor.Ycord - Source.Ycord
        d = np.sqrt(dx**2 + dy**2)
        capDistance = dx * self.SemiPerimeter / (m.pi * self.diameter) # É apenas uma aproximação, tem como calcular melhor esse valor

        case1 = (Source.Ycord + Sensor.Ycord + capDistance) <= d
        case2 = (2 * self.height - Source.Ycord + Sensor.Ycord + capDistance) <= d
        Cond1 = case1 or case2
        Cond2 = not (Source.OnCap or Sensor.OnCap)
        if Cond1 and Cond2:
            if Source.Ycord > SemiHeight:
                YAuxCord = self.height
            else:
                YAuxCord = 0

            AuxPoint1 = VesselPoint(0, YAuxCord, -2)
            AuxPoint2 = VesselPoint(0, YAuxCord, -2)
            self.__GenPlot = False  # Desabilitando temporariamente os gráficos para melhor desempenho

            def CalcVerticalClone(Xaux):
                """Função auxiliar a distância entre dois pontos no casco, mas de forma que o caminho passe pelo tampo.

                Arguments:
                Xaux {[float array]} -- [vetor com as coordenadas x dos pontos auxiliares]

                Returns:
                    [dist] -- [distância calculada]
                """

                AuxPoint1.SetXcord(Xaux[0])
                AuxPoint2.SetXcord(Xaux[1])
                self.__AuxCoords(AuxPoint1)
                self.__AuxCoords(AuxPoint2)
                AuxPoint1.SetOnCap(False)
                dist1 = self.__calcDist(Source, AuxPoint1)
                AuxPoint1.SetOnCap(True)
                AuxPoint2.SetOnCap(True)
                dist2 = self.__calcDist(AuxPoint1, AuxPoint2)
                AuxPoint2.SetOnCap(False)
                dist3 = self.__calcDist(AuxPoint2, Sensor)
                totalDist = dist1 + dist2 + dist3
                return totalDist

            SearchRange = (0, self.diameter * m.pi)
            InitGuess = opt.brute(
                CalcVerticalClone, (SearchRange, SearchRange), Ns=6)
            FinalSearch = opt.minimize(
                CalcVerticalClone, x0=InitGuess, method="BFGS")
            # print(FinalSearch) -- Resultado da minimização
            dist = FinalSearch.get("fun")
            MinPos1 = FinalSearch.get("x")[0]
            MinPos2 = FinalSearch.get("x")[1]
            AuxPoint1.SetXcord(MinPos1)
            AuxPoint2.SetXcord(MinPos2)

        else:
            dist = -1

        return dist

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord >= 0 and Xcord <= self.diameter * m.pi
        C2 = Ycord > - self.SemiPerimeter * \
            1.01 and Ycord < (self.height + self.SemiPerimeter) * 1.01
        if C1 and C2:
            self.__SensorID += 1
            ID = self.__SensorID
            SensorCoords = VesselPoint(Xcord, Ycord, ID)
            self.__AuxCoords(SensorCoords)
            self.SensorList.append(SensorCoords)
        else:
            print("As coordenadas deste ponto estão fora do vaso")

    def StructuredSensorDistribution(self, lines, sensorsInLine, x0, y0, dx, dy, aligned):
        for i in range(0, lines):
            if not aligned and (-1)**(i + 1) == 1:
                x1 = x0 + dx / 2
            else:
                x1 = x0
            y1 = y0 + i * dy
            for j in range(0, sensorsInLine):
                x = x1 + j * dx
                y = y1
                self.AddSensor(Xcord=x, Ycord=y)

    def FindFurthestPoint(self):
        # Revisar este função
        def CalcDistRemotePoint(x):
            distances = self.calcAllDist(SourceX=x[0], SourceY=x[1], IDs=[-1])
            return np.min(distances)

        def CallBack(xk):
            print(xk)
            maxDist = CalcDistRemotePoint(xk)
            print("Max distance: " + str(maxDist))

        BruteRes = opt.brute(lambda x: -CalcDistRemotePoint(x), ranges=[(
            0, self.diameter * m.pi), (-self.SemiPerimeter, self.height + self.SemiPerimeter)], Ns=5)
        # print(BruteRes)
        res = opt.minimize(lambda x: -CalcDistRemotePoint(x), method='L-BFGS-B', bounds=[(
            0, self.diameter * m.pi), (-self.SemiPerimeter, self.height + self.SemiPerimeter)], x0=BruteRes, callback=CallBack, options={'maxfun': 200, 'ftol': 0.0000001})

        return res.get('x')

    def __removeSensors(self, IDs):
        # Remove os sensores que não estão na lista de IDs
        if IDs == [-1]:
            pass
        else:
            ValidSensors = []
            InvalidSensors = []
            for sensor in self.SensorList:
                try:
                    IDs.index(sensor.ID)
                    ValidSensors.append(sensor)
                except:
                    InvalidSensors.append(sensor)

            self.SensorList = ValidSensors
            self.__tempSensorList = InvalidSensors

    def __returnSensors(self):
        # Os sensores sempre voltam ordenados às suas posições
        if not self.__tempSensorList == None:
            temp = self.SensorList + self.__tempSensorList
            self.SensorList = [None] * len(temp)
            for sensor in temp:
                self.SensorList[sensor.ID] = sensor

    def __initializeTimes(self):
        self.t_samecap = 0
        self.t_wallcap = 0
        self.t_wall = 0
        self.t_captocap = 0
        self.i_samecap = 0
        self.i_wallcap = 0
        self.i_wall = 0
        self.i_captocap = 0

    def __printTimes(self):
        print("\nTempos de cada modo de cálculo de distâncias")
        try:
            print("Mesmo tampo: " + str(round(self.t_samecap, 3)) +
                  " em " + str(self.i_samecap) + " avaliações / tempo médio: " + str(round(self.t_samecap / self.i_samecap, 4)) + " s")
        except:
            pass
        try:
            print("Corpo tampo: " + str(round(self.t_wallcap, 3)) +
                  " em " + str(self.i_wallcap) + " avaliações / tempo médio: " + str(round(self.t_wallcap / self.i_wallcap, 4)) + " s")
        except:
            pass
        try:
            print("Tampo tampo: " + str(round(self.t_captocap, 3)) +
                  " em " + str(self.i_captocap) + " avaliações / tempo médio: " + str(round(self.t_captocap / self.i_captocap, 4)) + " s")
        except:
            pass
        try:
            print("Corpo: " + str(round(self.t_wall, 3)) +
                  " em " + str(self.i_wall) + " avaliações / tempo médio: " + str(round(self.t_wall / self.i_wall, 4)) + " s")
        except:
            pass

        print("\n")

    def calcAllDist(self, SourceX, SourceY, IDs):
        """
        # Inicialização dos tempos acumulados de cálculo de distâncias  - medição de performance
        self.__initializeTimes()
        """

        self.__removeSensors(IDs)
        Source = VesselPoint(SourceX, SourceY, -1)
        self.__AuxCoords(Source)

        t0 = time.time()

        MinDistances = np.zeros(len(self.SensorList))
        i = -1
        for sensor in self.SensorList:
            i += 1

            distDirect = self.__calcDist(Source, sensor)
            distVClone = -1
            #distVClone = self.__DistVClone(Source, sensor)
            if distVClone == -1:
                distVClone = distDirect * 10
            else:
                print("Distância direta: " + str(distDirect))
                print("Clone vertical :" + str(distVClone))

            Distances = [distDirect, distVClone]
            MinDistances[i] = np.min(Distances)

        self.__returnSensors()

        t1 = time.time()
        dt = t1 - t0

        # Report dos tempos para calcular todas as distâncias
        """
        self.__printTimes()
        print("Tempo total: " + str(round(dt, 5)))
        """

        return MinDistances

    def __SimplifiedDistances(self, x, y, IDs):
        self.__removeSensors(IDs)
        distances = []
        for sensor in self.SensorList:
            dist1 = np.sqrt((x - sensor.Xcord)**2 + (y - sensor.Ycord)**2)
            dist2 = np.sqrt((x - sensor.Xcord + self.diameter *
                             m.pi)**2 + (y - sensor.Ycord)**2)
            dist3 = np.sqrt((x - sensor.Xcord - self.diameter *
                             m.pi)**2 + (y - sensor.Ycord)**2)
            minDist = np.min([dist1, dist2, dist3])
            distances.append(minDist)

        self.__returnSensors()
        return distances

    def returnDeltaT(self, x, y, IDs, mode):
        auxMode = self.CalcMode
        auxSection = self.SectionMode

        if mode == 'geodesic':
            self.setCalcMode('geodesic')
            self.__AllSensorsAuxCoords()
            distances = self.calcAllDist(x, y, IDs)
        elif mode == 'reg':
            self.setCalcMode('section')
            self.setSectionMode('reg')
            self.__AllSensorsAuxCoords()
            distances = self.calcAllDist(x, y, IDs)
        elif mode == 'inc':
            self.setCalcMode('section')
            self.setSectionMode('inc')
            self.__AllSensorsAuxCoords()
            distances = self.calcAllDist(x, y, IDs)
        elif mode == 'simple':
            "Modo simplificado - planificado"
            distances = self.__SimplifiedDistances(x, y, IDs)
        elif mode == 'original':
            "Usar modo atual"
            distances = self.calcAllDist(x, y, IDs)
        else:
            print("Modo inexistente")

        if mode != 'original' and mode != 'simple':
            self.setCalcMode(auxMode)
            self.setSectionMode(auxSection)
            self.__AllSensorsAuxCoords()

        NPdist = np.array(distances)
        times = (NPdist - NPdist[0]) / self.veloc
        
        return times

    def __orderMembers(self, TimesToSensors):
        type = [('ID', int), ('time', float)]
        data = np.array(TimesToSensors, dtype=type)
        OrderedMembers = np.sort(data, order='ID')

        return OrderedMembers

    def __orderByTime(self, TimesToSensors):
        dtype = [("ID", int), ("time", float)]
        NPtimes = np.array(TimesToSensors, dtype=dtype)
        NPtimes = np.sort(NPtimes, order="time")

        return NPtimes

    def __InitialKick(self, TimesToSensors):
        temp = self.__orderByTime(TimesToSensors)
        times = []
        Nsensors = 3
        for element in temp[:Nsensors]:
            (ID, time) = element
            times.append(time)

        weights = np.array(times / np.max(times))
        weights += 1

        x = 0
        y = 0
        j = 0
        sum = 0
        for element in temp[:Nsensors]:
            (ID, time) = element
            sensor = self.__getSensorbyID(ID)
            x += sensor.Xcord / weights[j]
            y += sensor.Ycord / weights[j]
            sum += 1 / weights[j]
            j += 1

        x0 = [x / sum, y / sum]
        """
        # Definindo como sendo o sensor mais próximo
        (ID, time) = temp[0]
        sensor = self.__getSensorbyID(ID)
        x0 = [sensor.Xcord, sensor.Ycord]
        """
        return (x0)

    def simpleLocation(self, TimesToSensors):
        x0 = self.__InitialKick(TimesToSensors)

        data = self.__orderMembers(TimesToSensors)
        IDs = []
        MeasTimes = []

        for member in data:
            (ID, time) = member
            IDs.append(ID)
            MeasTimes.append(time)

        MeasTimes = np.array(MeasTimes)
        normalizer = 1 / (np.max(MeasTimes) + 1)

        def CalcResidue(x):
            tcalc = self.returnDeltaT(x[0], x[1], IDs, 'simple')
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum(((tcalc - MeasTimes) * normalizer)**2))
            f = np.log10(residue)
            return f

        res = opt.minimize(CalcResidue, x0, method='BFGS')
        # print("Localização simplificada:")
        # print(res.get("x"))

        return res.get("x")

    def completeLocation(self, TimesToSensors):
        # Inicialização dos tempos acumulados
        #self.__initializeTimes()

        #x0 = self.__InitialKick(TimesToSensors)

        data = self.__orderMembers(TimesToSensors)
        (firstID, t0) = data[0]
        IDs = []
        MeasTimes = []

        for member in data:
            (ID, TOF) = member
            IDs.append(ID)
            MeasTimes.append(TOF - t0)

        MeasTimes = np.array(MeasTimes)
        normalizer = 1 / (np.max(MeasTimes) + 1)

        def CalcResidue(x):
            tcalc = self.returnDeltaT(x[0], x[1], IDs, 'original')
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum(((tcalc - MeasTimes) * normalizer)**2))
            f = np.log10(residue)

            return f

        # options={"gtol": 1E-4}
        bounds = [(-0.01 * self.diameter * m.pi, 1.01 * self.diameter * m.pi),
                  (-1.01 * self.SemiPerimeter, 1.01 * (self.height + self.SemiPerimeter))]
        maxiter = 1000
        polish = False

        """
        res = opt.minimize(CalcResidue, x0=x0, method='L-BFGS-B', options={"gtol": 1E-4}, bounds=bounds)
        """

        res = opt.differential_evolution(
            CalcResidue, bounds=bounds, maxiter=maxiter, polish=polish)


        print(res)  # - Resultado da otimização

        #self.__printTimes()

        return res.get("x")

    def fCostMap(self, TimesToSensors):
        data = self.__orderMembers(TimesToSensors)
        IDs = []
        MeasTimes = []

        for member in data:
            (ID, time) = member
            IDs.append(ID)
            MeasTimes.append(time)

        MeasTimes = np.array(MeasTimes)

        def CalcResidue(x):
            tcalc = self.returnDeltaT(x[0], x[1], IDs, 'original')
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum((tcalc - MeasTimes)**2))
            return residue

        N = 20
        """
        x_array = np.linspace(0, self.diameter * m.pi, num=N)
        y_array = np.linspace(-self.SemiPerimeter,
                              self.height + self.SemiPerimeter, num=N)
        """

        x_array = np.linspace(0, self.diameter * m.pi, num=N)
        y_array = np.linspace(self.height * 0, self.height +
                              self.SemiPerimeter * 0, num=N)

        map = np.zeros((N, N))
        i = 0
        for x in x_array:
            j = 0
            for y in y_array:
                map[i, j] = CalcResidue([x, y])
                print("i: " + str(i) + " - j: " + str(j))
                j += 1
            i += 1

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        x_array, y_array = np.meshgrid(x_array, y_array)

        surf = ax.plot_surface(x_array, y_array, map, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        plt.show()
