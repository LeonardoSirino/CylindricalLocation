import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math as m
import numpy as np


class VesselPoint():
    """Classe para definir as propriedades de um ponto no vaso
    """

    def __init__(self, Xcord, Ycord, ID):
        self.Xcord = Xcord
        self.Ycord = Ycord
        self.ID = ID
        self.Valid = True

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

    def SetValid(self, Valid):
        """Validade de um ponto para se calcular a distância.
        Os clones horizontais de pontos no tampo não são válidos e portanto não se deve a calcular a ditância até eles.
        Esses são mantidos para manter a equivalência entre os vetores de sensores e clones

        Arguments:
            Valid {[boolean]} -- [validade do ponto]
        """

        self.Valid = Valid

    def __str__(self):
        if self.Valid:
            Valido = "Valido"
        else:
            Valido = "Invalido"
        text = Valido + "ID: " + \
            str(self.ID) + " x: " + str(self.Xcord) + " y: " + str(self.Ycord)
        return text


class CylindricalLocation():
    def __init__(self, diameter, height):
        self.diameter = diameter
        self.height = height
        self.SensorList = []
        self.veloc = 1
        self.__SensorListHClone = []
        self.__Xpath = []
        self.__Ypath = []
        self.__GenPlot = False
        self.__PointsonPlot = 100
        self.__SensorID = -1
        self.__tempSensorList = None
        self.__tempSensorListHClone = None
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

    def setCalcMode(self, mode):
        self.CalcMode = mode

    def setSectionMode(self, mode):
        self.SectionMode = mode

    def set_f(self, f):
        self.f = f
        self.cap = geo.Geodesic(self.diameter / 2, f)
        result = self.cap.Inverse(lat1=0, lon1=0, lat2=90, lon2=0)
        self.SemiPerimeter = result.get("s12")
        self.__DrawVessel()

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
        self.__DrawVessel()

    def SetVelocity(self, velocity):
        "Definição da velocidade em mm/s"
        self.veloc = velocity

    def PrintAllSensors(self):
        print("Sensores originais:")
        for sensor in self.SensorList:
            print(sensor)

        print("\n Clones horizontais:")
        for sensor in self.__SensorListHClone:
            print(sensor)

    def __getSensorbyID(self, ID):
        for sensor in self.SensorList:
            if sensor.ID == ID:
                return sensor

        return None

    def __DrawVessel(self):
        self.__VesselXPath = self.__PlotRectangle(0, 0, self.diameter * m.pi, -self.SemiPerimeter).get("xpath") + self.__PlotRectangle(
            0, 0, self.diameter * m.pi, self.height).get("xpath") + self.__PlotRectangle(0, self.height, self.diameter * m.pi, self.SemiPerimeter).get("xpath")
        self.__VesselYPath = self.__PlotRectangle(0, 0, self.diameter * m.pi, -self.SemiPerimeter).get("ypath") + self.__PlotRectangle(
            0, 0, self.diameter * m.pi, self.height).get("ypath") + self.__PlotRectangle(0, self.height, self.diameter * m.pi, self.SemiPerimeter).get("ypath")

    def __PlotRectangle(self, x0, y0, width, height):
        xpath = [x0, x0 + width, x0 + width, x0, x0]
        ypath = [y0, y0, y0 + height, y0 + height, y0]
        return {"xpath": xpath, "ypath": ypath}

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

    def __PlotonCap(self, lat1, lat2, lon1, lon2, cap):
        if self.CalcMode == 'geodesic':
            smax = self.cap.Inverse(lat1=lat1, lat2=lat2,
                                    lon1=lon1, lon2=lon2).get("s12")
            Path = self.cap.InverseLine(
                lat1=lat1, lat2=lat2, lon1=lon1, lon2=lon2)
            lenghts = np.linspace(0, smax, self.__PointsonPlot)
            oldpathDiam = 0
            for s in lenghts:
                position = Path.Position(s12=s)
                pathDiam = (position.get("lon2") + 180) * \
                    self.diameter * m.pi / 360
                if abs(oldpathDiam - pathDiam) > self.diameter * m.pi * 0.8:
                    self.__Ypath.append(m.nan)
                    self.__Xpath.append(m.nan)
                oldpathDiam = pathDiam
                lataux = position.get("lat2")
                Saux = self.cap.Inverse(
                    lat1=0, lon1=0, lat2=lataux, lon2=0).get("s12")
                if cap == "sup":
                    Saux = Saux + self.height
                else:
                    Saux = -Saux
                self.__Ypath.append(Saux)
                self.__Xpath.append(pathDiam)
        elif self.CalcMode == 'section':
            print("Não há plots no tampo para esse modo")
        else:
            print("Modo inválido")

    def __PlotonWall(self, x1, x2, y1, y2):
        xpath = np.linspace(x1, x2, self.__PointsonPlot)
        xpath = xpath.tolist()
        ypath = np.linspace(y1, y2, self.__PointsonPlot)
        ypath = ypath.tolist()
        # Correção da coordenada X
        aux = True
        xcorrect = []
        ycorrect = []
        k = 0
        for x in xpath:
            if x > self.diameter * m.pi:
                x = x - self.diameter * m.pi
                if aux:
                    xcorrect.append(m.nan)
                    ycorrect.append(m.nan)
                    aux = False
            elif x < 0:
                x = x + self.diameter * m.pi
                if aux:
                    xcorrect.append(m.nan)
                    ycorrect.append(m.nan)
                    aux = False
            xcorrect.append(x)
            ycorrect.append(ypath[k])
            k += 1
        self.__Xpath += xcorrect
        self.__Ypath += ycorrect

    def __calcDist(self, P1, P2):
        """
        Função para calcular distância entre dois pontos em posições quaisquer do vaso
        """
        if P1.OnCap and P2.OnCap:
            if P1.Cap == P2.Cap:  # Dois pontos no mesmo tampo
                dist = self.__DistSameCap(P1, P2)
            else:  # Pontos em tampos opostos
                if P1.Cap == "sup":
                    Psup = P1
                    Pinf = P2
                else:
                    Psup = P2
                    Pinf = P1
                dist = self.__DistCaptoCap(Psup, Pinf)
        # Distância entre um ponto no casco e outro no tampo
        elif P1.OnCap ^ P2.OnCap:
            dist = self.__DistWalltoCap(P1, P2)
        else:  # Distãncia entre pontos no casco
            dist = m.sqrt((P1.Xcord - P2.Xcord) **
                          2 + (P1.Ycord - P2.Ycord)**2)
            if self.__GenPlot:
                self.__PlotonWall(x1=P1.Xcord, x2=P2.Xcord,
                                  y1=P1.Ycord, y2=P2.Ycord)

        return dist

    def __DistSameCap(self, P1, P2):
        if self.CalcMode == 'geodesic':
            res = self.cap.Inverse(
                lat1=P1.Lat, lat2=P2.Lat, lon1=P1.Lon, lon2=P2.Lon)
            dist = res.get("s12")
            if self.__GenPlot:
                self.__PlotonCap(lat1=P1.Lat, lat2=P2.Lat,
                                 lon1=P1.Lon, lon2=P2.Lon, cap=P1.Cap)
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

                theta1 = np.arccos(np.dot(v1c, v12) /
                                   (np.linalg.norm(v1c) * np.linalg.norm(v12)))

                theta2 = np.arccos(np.dot(v2c, v12) /
                                   (np.linalg.norm(v2c) * np.linalg.norm(v12)))

                if theta1 > m.pi / 2:
                    u1 = -u1

                if theta2 > m.pi / 2:
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
        AuxGenPlot = self.__GenPlot
        self.__GenPlot = False  # Desabilitando temporariamente os gráficos para melhor desempenho

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
            dist1 = self.__calcDist(Pwall, AuxPoint)
            AuxPoint.SetOnCap(True)
            dist2 = self.__calcDist(AuxPoint, Pcap)
            totalDist = dist1 + dist2
            return totalDist

        InitGuess = opt.brute(
            CalcWallToCap, ((0, m.pi * self.diameter),), Ns=6)
        FinalSearch = opt.minimize(CalcWallToCap, x0=InitGuess, method="BFGS")
        # print(FinalSearch) -- Resultado da minimização
        dist = FinalSearch.get("fun")
        MinPos = FinalSearch.get("x")[0]
        AuxPoint.SetXcord(MinPos)
        BestPoint = AuxPoint

        self.__GenPlot = AuxGenPlot

        if self.__GenPlot:  # Plot do caminho que passa pelo casco e pelo tampo
            self.__AuxCoords(BestPoint)
            BestPoint.SetOnCap(False)
            # Plot do caminho no casco
            dist1 = self.__calcDist(Pwall, BestPoint)
            BestPoint.SetOnCap(True)
            # Plot do caminho no tampo
            dist2 = self.__calcDist(BestPoint, Pcap)

        return dist

    def __DistCaptoCap(self, Psup, Pinf):
        AuxPoint1 = VesselPoint(0, self.height, -2)
        AuxPoint2 = VesselPoint(0, 0, -2)
        AuxGenPlot = self.__GenPlot
        self.__GenPlot = False  # Desabilitando temporariamente os gráficos para melhor desempenho

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

        # -- Chute inicial com otimização bruta
        SearchRange = (0, m.pi * self.diameter)
        InitGuess = opt.brute(CalcCaptoCap, (SearchRange, SearchRange), Ns=6)
        FinalSearch = opt.minimize(CalcCaptoCap, x0=InitGuess, method="BFGS")
        # print(FinalSearch)  # -- Resultado da minimização
        dist = FinalSearch.get("fun")
        MinPosSup = FinalSearch.get("x")[0]
        MinPosInf = FinalSearch.get("x")[1]
        AuxPoint1.SetXcord(MinPosSup)
        AuxPoint2.SetXcord(MinPosInf)

        self.__GenPlot = AuxGenPlot

        if self.__GenPlot:  # Plot do caminho que passa pelo casco e pelo tampo
            self.__AuxCoords(AuxPoint1)
            AuxPoint1.SetOnCap(True)
            # Plot do caminho no tampo superior
            dist1 = self.__calcDist(Psup, AuxPoint1)
            AuxPoint1.SetOnCap(False)
            AuxPoint2.SetOnCap(False)
            # Plot do caminho no casco
            dist2 = self.__calcDist(AuxPoint1, AuxPoint2)
            AuxPoint2.SetOnCap(True)
            # Plot do caminho no tampo inferior
            dist3 = self.__calcDist(AuxPoint2, Pinf)

        return dist

    def __DistVClone(self, Source, Sensor):
        SemiHeight = self.height / 2
        Cond1 = (Source.Ycord > SemiHeight and Sensor.Ycord > SemiHeight) or (
            Source.Ycord < SemiHeight and Sensor.Ycord < SemiHeight)
        Cond2 = not (Source.OnCap or Sensor.OnCap)
        if Cond1 and Cond2:
            if Source.Ycord > SemiHeight:
                YAuxCord = self.height
            else:
                YAuxCord = 0

            AuxPoint1 = VesselPoint(0, YAuxCord, -2)
            AuxPoint2 = VesselPoint(0, YAuxCord, -2)
            AuxGenPlot = self.__GenPlot
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

            self.__GenPlot = AuxGenPlot

            if self.__GenPlot:  # Plot do caminho que passa pelo casco e pelo tampo
                self.__AuxCoords(AuxPoint1)
                AuxPoint1.SetOnCap(False)
                # Plot do caminho da fonte até o ponto auxiliar 1
                dist1 = self.__calcDist(Source, AuxPoint1)
                AuxPoint1.SetOnCap(True)
                AuxPoint2.SetOnCap(True)
                # Plot do caminho entre pontos auxiliares
                dist2 = self.__calcDist(AuxPoint1, AuxPoint2)
                AuxPoint2.SetOnCap(False)
                # Plot do ponto auxiliar 2 até o sensor
                dist3 = self.__calcDist(AuxPoint2, Sensor)

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
            SensorHCloneCoords = VesselPoint(0, Ycord, ID)
            if SensorCoords.Ycord >= 0 and SensorCoords.Ycord <= self.height:  # Ponto no casco do vaso
                # Sensor posicionado na porção direita do casco
                if SensorCoords.Xcord >= self.diameter * m.pi / 2:
                    XHClone = SensorCoords.Xcord - self.diameter * m.pi
                else:
                    # Sensor posicionado na porção esquerda do casco
                    XHClone = SensorCoords.Xcord + self.diameter * m.pi

                SensorHCloneCoords.SetXcord(XHClone)
            else:  # Ponto em um dos tampos do vaso
                SensorHCloneCoords.SetXcord(m.nan)
                SensorHCloneCoords.SetValid(False)

            self.__AuxCoords(SensorHCloneCoords)
            self.__SensorListHClone.append(SensorHCloneCoords)
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
        def CalcDistRemotePoint(x):
            distances = self.calcAllDist(
                SourceX=x[0], SourceY=x[1], GenPlot=False)
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

            ValidSensors = []
            InvalidSensors = []
            for sensor in self.__SensorListHClone:
                try:
                    IDs.index(sensor.ID)
                    ValidSensors.append(sensor)
                except:
                    InvalidSensors.append(sensor)

            self.__SensorListHClone = ValidSensors
            self.__tempSensorListHClone = InvalidSensors

    def __returnSensors(self):
        # Os sensores sempre voltam ordenados às suas posições
        if not self.__tempSensorList == None:
            temp = self.SensorList + self.__tempSensorList
            self.SensorList = [None] * len(temp)
            for sensor in temp:
                self.SensorList[sensor.ID] = sensor

            temp = self.__SensorListHClone + self.__tempSensorListHClone
            self.__SensorListHClone = [None] * len(temp)
            for sensor in temp:
                self.__SensorListHClone[sensor.ID] = sensor

    def calcAllDist(self, SourceX, SourceY, GenPlot, IDs):
        self.__removeSensors(IDs)
        self.__GenPlot = GenPlot
        Source = VesselPoint(SourceX, SourceY, -1)
        self.__AuxCoords(Source)
        MinDistances = []
        if GenPlot:
            plt.plot(self.__VesselXPath, self.__VesselYPath)
            SensorX = []
            SensorY = []
            for sensor in self.SensorList:
                SensorX.append(sensor.Xcord)
                SensorY.append(sensor.Ycord)
            SensorX.append(SourceX)
            SensorY.append(SourceY)
            plt.plot(SensorX, SensorY, ".")

        i = -1
        for sensor in self.SensorList:
            i += 1
            # Limpando o histórico do plot para cada sensor
            self.__Xpath = []
            self.__Ypath = []

            distDirect = self.__calcDist(Source, sensor)
            self.__Xpath.append(-10)
            self.__Ypath.append(m.nan)

            HClloneSensor = self.__SensorListHClone[i]
            if HClloneSensor.Valid:
                distHClone = self.__calcDist(Source, HClloneSensor)
            else:
                distHClone = distDirect * 10
                self.__Xpath.append(0)
                self.__Ypath.append(0)
            self.__Xpath.append(-10)
            self.__Ypath.append(m.nan)

            distVClone = self.__DistVClone(Source, sensor)
            if distVClone == -1:
                distVClone = distDirect * 10
                self.__Xpath.append(0)
                self.__Ypath.append(0)
            self.__Xpath.append(-10)
            self.__Ypath.append(m.nan)

            Distances = [distDirect, distHClone, distVClone]
            MinDistances.append(np.min(Distances))

            if GenPlot:
                MinDistIndex = Distances.index(np.min(Distances))
                XAllPath = np.array(self.__Xpath)
                YAllPath = np.array(self.__Ypath)
                Xpaths = np.split(XAllPath, np.where(XAllPath == -10)[0])
                Ypaths = np.split(YAllPath, np.where(XAllPath == -10)[0])
                XBestPath = Xpaths[MinDistIndex]
                XBestPath[0] = m.nan
                YBestPath = Ypaths[MinDistIndex]
                YBestPath[0] = m.nan
                if MinDistIndex == 0:
                    path = "CD"
                elif MinDistIndex == 1:
                    path = "CH"
                else:
                    path = "CV"

                print("Fonte: " + str(round(Source.Xcord, 1)) + " - " + str(round(Source.Ycord, 1)) + " -- Sensor: " + str(round(
                    sensor.Xcord, 1)) + " - " + str(round(sensor.Ycord, 1)) + " -- Distância: " + str(round(np.min(Distances), 3)) + " -- " + path)

                plt.plot(XBestPath, YBestPath)

        if GenPlot:
            plt.show()

        self.__returnSensors()
        return MinDistances

    def __SimplifiedDistances(self, x, y, IDs):
        self.__removeSensors(IDs)
        distances = []
        for (sensor, clone) in zip(self.SensorList, self.__SensorListHClone):
            dist1 = np.sqrt((x - sensor.Xcord)**2 + (y - sensor.Ycord)**2)
            if clone.Valid:
                dist2 = np.sqrt((x - clone.Xcord)**2 + (y - clone.Ycord)**2)
            else:
                dist2 = 10 * dist1
            distances.append(np.min([dist1, dist2]))

        self.__returnSensors()
        return distances

    def returnDeltaT(self, x, y, IDs, exact):
        if exact:
            distances = self.calcAllDist(x, y, False, IDs)
        else:
            distances = self.__SimplifiedDistances(x, y, IDs)
        NPdist = np.array(distances)
        NPdist += -np.min(NPdist)
        times = NPdist / self.veloc
        return times

    def __orderMembers(self, TimesToSensors):
        IDs = []
        for member in TimesToSensors:
            (ID, time) = member
            IDs.append(ID)

        IDs.sort()
        OrderedMembers = [0] * len(IDs)
        for member in TimesToSensors:
            (ID, time) = member
            i = IDs.index(ID)
            OrderedMembers[i] = member

        return OrderedMembers

    def __orderByTime(self, TimesToSensors):
        dtype = [("ID", int), ("time", float)]
        NPtimes = np.array(TimesToSensors, dtype=dtype)
        NPtimes = np.sort(NPtimes, order="time")

        return NPtimes

    def simpleLocation(self, TimesToSensors):
        temp = self.__orderByTime(TimesToSensors)
        times = []
        for element in temp[:3]:
            (ID, time) = element
            times.append(time)

        weights = np.array(times / np.max(times))
        weights += 1

        x = 0
        y = 0
        j = 0
        sum = 0
        for element in temp[:3]:
            (ID, time) = element
            sensor = self.__getSensorbyID(ID)
            x += sensor.Xcord / weights[j]
            y += sensor.Ycord / weights[j]
            sum += 1 / weights[j]
            j += 1
            
        x0 = [x / sum, y / sum]
        print("x0")
        print(x0)

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
            tcalc = self.returnDeltaT(x[0], x[1], IDs, False)
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum(((tcalc - MeasTimes) * normalizer)**2))
            return residue

        res = opt.minimize(CalcResidue, x0, method='BFGS')
        # print("Localização simplificada:")
        # print(res.get("x"))

        return res.get("x")

    def completeLocation(self, TimesToSensors):
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
            tcalc = self.returnDeltaT(x[0], x[1], IDs, True)
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum(((tcalc - MeasTimes) * normalizer)**2))
            return residue

        x0 = self.simpleLocation(TimesToSensors)

        res = opt.minimize(CalcResidue, x0=x0, method='L-BFGS-B', options={"gtol": 3E-3}, bounds=[
                           (-0.01 * self.diameter * m.pi, 1.01 * self.diameter * m.pi), (-1.01 * self.SemiPerimeter, 1.01 * (self.height + self.SemiPerimeter))])
        """
        res = opt.differential_evolution(CalcResidue, bounds=[
                                         (-0.01 * self.diameter * m.pi, 1.01 * self.diameter * m.pi), (-1.01 * self.SemiPerimeter, 1.01 * (self.height + self.SemiPerimeter))])
        """

        print(res)  # - Resultado da otimização

        return res.get("x")
