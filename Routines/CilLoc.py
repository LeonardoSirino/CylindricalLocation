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
        for s in lenghts:
            position = Path.Position(s12=s)
            pathDiam = (position.get("lon2") + 180) * \
                self.diameter * m.pi / 360
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
                # Definir os pontos corretos!
                dist = self.__DistCaptoCap(P1, P2)
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
        print("Pontos em tampos opostos")
        dist = 0
        return dist

    def __DistVClone(self, Source, Sensor):  # PlaceHolder
        dist = 0
        return dist

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord > 0 and Xcord < self.diameter * m.pi
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
        if GenPlot:
            plt.plot(self.__VesselXPath, self.__VesselYPath)
            SensorX = []
            SensorY = []
            for sensor in self.SensorList:
                SensorX.append(sensor.Xcord)
                SensorY.append(sensor.Ycord)
            for sensor in self.__SensorListHClone:
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
            dist = self.__calcDist(Source, sensor)
            print(str(i) + " - Distância direta: " + str(dist))
            self.__Xpath.append(m.nan)
            self.__Ypath.append(m.nan)
            HClloneSensor = self.__SensorListHClone[i]
            if HClloneSensor.Valid:
                dist = self.__calcDist(Source, HClloneSensor)
                print(str(i) + " - Distância clone horizontal: " + str(dist))
                self.__Xpath.append(m.nan)
                self.__Ypath.append(m.nan)
            if GenPlot:
                plt.plot(self.__Xpath, self.__Ypath)
            i += 1

        if GenPlot:
            plt.show()
