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

    def set_f(self, f):
        self.f = f
        self.cap = geo.Geodesic(self.diameter / 2, f)
        result = self.cap.Inverse(lat1=0, lon1=0, lat2=90, lon2=0)
        self.SemiPerimeter = result.get("s12")
        self.__DrawVessel()

    def set_semiPerimeter(self, SemiPerimeter):
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

    def __DrawVessel(self):
        self.__VesselXPath = self.__PlotRectangle(0, 0, self.diameter * m.pi, -self.SemiPerimeter).get("xpath") + self.__PlotRectangle(
            0, 0, self.diameter * m.pi, self.height).get("xpath") + self.__PlotRectangle(0, self.height, self.diameter * m.pi, self.SemiPerimeter).get("xpath")
        self.__VesselYPath = self.__PlotRectangle(0, 0, self.diameter * m.pi, -self.SemiPerimeter).get("ypath") + self.__PlotRectangle(
            0, 0, self.diameter * m.pi, self.height).get("ypath") + self.__PlotRectangle(0, self.height, self.diameter * m.pi, self.SemiPerimeter).get("ypath")

    def __PlotRectangle(self, x0, y0, width, height):
        xpath = [x0, x0 + width, x0 + width, x0, x0]
        ypath = [y0, y0, y0 + height, y0 + height, y0]
        return {"xpath": xpath, "ypath": ypath}

    def __AuxCoords(self, Coords):
        """
        Função para calcular as coordenadas auxliares (latitude e longitude) quando a coordenada estiver no tampo
        """
        if Coords.Ycord >= self.height or Coords.Ycord <= 0:
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
            Coords.SetOnCap(True)
        else:
            Coords.SetOnCap(False)

    def __PlotonCap(self, lat1, lat2, lon1, lon2, cap):
        smax = self.cap.Inverse(lat1=lat1, lat2=lat2,
                                lon1=lon1, lon2=lon2).get("s12")
        Path = self.cap.InverseLine(lat1=lat1, lat2=lat2, lon1=lon1, lon2=lon2)
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
                res = self.cap.Inverse(
                    lat1=P1.Lat, lat2=P2.Lat, lon1=P1.Lon, lon2=P2.Lon)
                dist = res.get("s12")
                if self.__GenPlot:
                    self.__PlotonCap(lat1=P1.Lat, lat2=P2.Lat,
                                     lon1=P1.Lon, lon2=P2.Lon, cap=P1.Cap)

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

    def simpleLocation(self, TimesToSensors):
        """
        Melhorar Chute inicial e definir pesos nos resíduos
        Passar resultado desse processo para a localização completa (considerando as geodésicas)
        """

        data = self.__orderMembers(TimesToSensors)
        IDs = []
        MeasTimes = []

        for member in data:
            (ID, time) = member
            IDs.append(ID)
            MeasTimes.append(time)

        MeasTimes = np.array(MeasTimes)
        normalizer = 1 / np.max(MeasTimes)

        def CalcResidue(x):
            tcalc = self.returnDeltaT(x[0], x[1], IDs, False)
            tcalc = np.array(tcalc)
            residue = np.sqrt(np.sum(((tcalc - MeasTimes) * normalizer)**2))
            return residue

        res = opt.minimize(CalcResidue, x0=[400, 500], method='BFGS')
        print(res)
        """
        x = res.get('x')
        print(MeasTimes)
        tcalc = self.returnDeltaT(x[0], x[1], IDs, False)
        print(tcalc)
        print((MeasTimes - tcalc)**2)
        residue = CalcResidue(x)
        print(residue)
        """
