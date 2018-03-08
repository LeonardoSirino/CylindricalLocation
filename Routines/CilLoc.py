import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math as m
import numpy as np


class VesselPoint():
    """Classe para definir as propriedades de um ponto no vaso
    """

    def __init__(self, Xcord, Ycord):
        self.Xcord = Xcord
        self.Ycord = Ycord
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


class CylindricalLocation():
    def __init__(self, diameter, height, SemiPerimeter):
        self.diameter = diameter
        self.height = height
        self.SemiPerimeter = SemiPerimeter
        self.f = SemiPerimeter * 4 * 2 / \
            (m.pi * diameter) - 1  # Achatamento do tampo
        self.SensorList = []
        self.__SensorListHClone = []
        self.cap = geo.Geodesic(self.diameter / 2, self.f)
        self.__Xpath = []
        self.__Ypath = []
        self.__GenPlot = False
        self.__PointsonPlot = 100
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

        minX = np.min([Pcap.Xcord, Pwall.Xcord])
        maxX = np.max([Pcap.Xcord, Pwall.Xcord])
        AuxPoint = VesselPoint(0, Yaux)
        SearchRange = slice(minX, maxX, (maxX - minX) / 100)
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

        InitGuess = opt.brute(CalcWallToCap, (SearchRange,))
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

    def __DistCaptoCap(self, Psup, Pinf):  # PlaceHolder
        minX = 0
        maxX = self.diameter * m.pi
        AuxPoint1 = VesselPoint(0, self.height)
        AuxPoint2 = VesselPoint(0, 0)
        SearchRange = slice(minX, maxX, (maxX - minX) / 10)
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

        InitGuess = opt.brute(CalcCaptoCap, (SearchRange, SearchRange))
        FinalSearch = opt.minimize(CalcCaptoCap, x0=InitGuess, method="BFGS")
        # print(FinalSearch) -- Resultado da minimização
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

    def __DistVClone(self, Source, Sensor):  # PlaceHolder
        SemiHeight = self.height/2
        Cond1 = (Source.Ycord > SemiHeight and Sensor.Ycord > SemiHeight) or (Source.Ycord < SemiHeight and Sensor.Ycord < SemiHeight)
        Cond2 = not (Source.OnCap or Sensor.OnCap)
        if Cond1 and Cond2:
            if Source.Ycord > SemiHeight:
                YAuxCord = self.height
            else:
                YAuxCord = 0

            minX = 0
            maxX = self.diameter * m.pi
            AuxPoint1 = VesselPoint(0, YAuxCord)
            AuxPoint2 = VesselPoint(0, YAuxCord)
            SearchRange = slice(minX, maxX, (maxX - minX) / 10)
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

            InitGuess = opt.brute(CalcVerticalClone, (SearchRange, SearchRange))
            FinalSearch = opt.minimize(CalcVerticalClone, x0=InitGuess, method="BFGS")
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
            dist=-1

        return dist

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord >= 0 and Xcord <= self.diameter * m.pi
        C2 = Ycord > - self.SemiPerimeter * \
            1.01 and Ycord < (self.height + self.SemiPerimeter) * 1.01
        if C1 and C2:
            SensorCoords = VesselPoint(Xcord, Ycord)
            self.__AuxCoords(SensorCoords)
            self.SensorList.append(SensorCoords)
            SensorHCloneCoords = VesselPoint(0, Ycord)
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

    def calcAllDist(self, SourceX, SourceY, GenPlot):
        self.__GenPlot = GenPlot
        Source = VesselPoint(SourceX, SourceY)
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


        i = 0
        for sensor in self.SensorList:
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
                print("Fonte: " + str(round(Source.Xcord,1)) + " - " + str(round(Source.Ycord,1)) + " -- Sensor: " + str(round(sensor.Xcord,1)) + " - " + str(round(sensor.Ycord,1)) + " -- Distância: " + str(round(np.min(Distances),3)))

                plt.plot(XBestPath, YBestPath)

        if GenPlot:
            plt.show()

        return MinDistances
